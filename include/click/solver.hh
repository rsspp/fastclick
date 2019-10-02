#ifndef SOLVER_H
#define SOLVER_H
#include <float.h>
#include <click/vector.hh>
#include <queue>
#include <vector>
#include <limits.h>
#include <algorithm>

inline double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}

struct cref {
    float load;
    int id;
};
class Compare
{
public:
    bool operator() (cref a, cref b)
    {
        return a.load < b.load;
    }
};

/**
 * Take some buckets of a list to make the load of a set of CPU reach target
 */
class BucketMapTargetProblem
{ public:
    Vector<int> transfer; //existing core id for each buckets
    Vector<float> target; //Target for destination
    Vector<float> max; //Max load to take from overloaded
    Vector<int> buckets_max_idx; //Load for each bucket
    Vector<float> buckets_load; //Load for each bucket
    float square_imbalance = 0; //Final load imbalance (after calling solve)

    BucketMapTargetProblem(int nbuckets, int nunderloaded, int noverloaded ) : min_cost(FLT_MAX) {
        transfer.resize(nbuckets,-1);
        buckets_load.resize(nbuckets);
        buckets_max_idx.resize(nbuckets);
        target.resize(nunderloaded);
        max.resize(noverloaded);
    }

    void solve(DeviceBalancer* balancer) {
        click_chatter("Assigning %d buckets to %d targets :", buckets_load.size(), target.size());

        if (unlikely(balancer->_verbose > 1)) {
            for (int i = 0; i < max.size(); i++) {
                click_chatter("Overloaded %d , %f",i, max[i]);
            }
            for (int i = 0; i < target.size(); i++) {
                click_chatter("Underloaded %d , %f",i, target[i]);
            }
            for (int i = 0; i < buckets_load.size(); i++) {
                click_chatter("Bucket %d , %f, oid %d",i, buckets_load[i], buckets_max_idx[i]);
            }
        }

        float target_imbalance = 0.01;//(o_load - u_loadà ;
        //Build priority for each overloaded
        //

        auto cmplt = [](cref left, cref right) { return left.load > right.load; };
        auto cmpgt = [](cref left, cref right) { return left.load < right.load; };

        float overload_allowed =  -0.01; //Take just not enough load of overloaded
        float underload_allowed = -0.01; //Give just not enough to underloaded
        float imbalance_u = 0;
        float imbalance_o = 0;
        float oa_min;
        float ua_min;
        float bottom_oa = overload_allowed;
        float top_oa;
        float bottom_ua = underload_allowed;
        float top_ua;

        float last_sq = FLT_MAX;
        float min_sq = FLT_MAX;
        float bottom_sq = FLT_MAX;
        float top_sq = FLT_MAX;
        int run = 1;
        float m = -1;
        int phase;
#define max_runs 10

        while(true) {

            Timestamp run_begin = Timestamp::now_steady();
            imbalance_u = 0;
            imbalance_o = 0;
            square_imbalance = 0;
            typedef std::priority_queue<cref, std::vector<cref>, Compare> bstack;
            std::vector<bstack> bstacks;
            bstacks.resize(max.size(),bstack());
            for (int i = 0; i < buckets_load.size(); i++) {
                bstacks[buckets_max_idx[i]].push(cref{buckets_load[i],i});
                transfer[i] = -1;
            }

            std::priority_queue<cref, std::vector<cref>, decltype(cmpgt)> overloaded(cmpgt);
            for (int i = 0; i < max.size(); i++) {
                overloaded.push(cref{max[i],i});
            }

            std::priority_queue<cref, std::vector<cref>, decltype(cmpgt)> underloaded(cmpgt);
            for (int i = 0; i < target.size(); i++) {
                underloaded.push(cref{target[i],i});
            }


            while (!overloaded.empty() && !underloaded.empty()) {
                next_core:
                //Select most overloaded core
                cref o = overloaded.top();
                overloaded.pop(); //Will be added bacj

                //Select biggest bucket
                bstack& buckets = bstacks[o.id];
                cref bucket = buckets.top();
                buckets.pop();//Bucket is removed forever

                //Assign its most used buckets
                std::vector<cref> save;
                while (!underloaded.empty()) {
                    cref u = underloaded.top();
                    underloaded.pop();

                    if (unlikely(balancer->_verbose > 2))
                        click_chatter("U%d load %f",u.id, u.load);
                    if (bucket.load < u.load + underload_allowed) {
                        u.load -= bucket.load;
                        o.load -= bucket.load;
                        transfer[bucket.id] = u.id;

                        if (unlikely(balancer->_verbose > 1))
                            click_chatter("Bucket %d to ucore %d, bucket load %f", bucket.id, u.id, bucket.load);
                        if (u.load > target_imbalance) {
                            underloaded.push(u);
                        } else {
                            imbalance_u += abs(u.load);
                            square_imbalance += u.load * u.load;
                            if (unlikely(balancer->_verbose > 1))
                                    click_chatter("Underloaded core %d is now okay with %f load", u.id, u.load);
                        }
                        goto bucket_assigned;
                    } else {
                        save.push_back(u);
                    }
                }

                if (unlikely(balancer->_verbose > 2))
                    click_chatter("Bucket %d UNMOVED, load %f", bucket.id, o.id, bucket.load);

                bucket_assigned:
                while (!save.empty()) {
                    underloaded.push(save.back());
                    save.pop_back();
                }

                if (o.load > - overload_allowed
                        && !buckets.empty()) {
                    overloaded.push(o);
                } else {
                    imbalance_o += abs(o.load);
                    square_imbalance += o.load * o.load;
                    if (unlikely(balancer->_verbose > 1))
                        click_chatter("Overloaded core %d is now okay with %f load. Empty : %d", o.id, o.load,buckets.empty());
                }
            }
            while (!overloaded.empty()) {
                auto o = overloaded.top();
                imbalance_o += abs(o.load);
                square_imbalance += o.load * o.load;
                overloaded.pop();
            }
            while (!underloaded.empty()) {
                auto u = underloaded.top();
                imbalance_u += abs(u.load);
                square_imbalance += u.load * u.load;
                underloaded.pop();
            }


            Timestamp run_end = Timestamp::now_steady();
            unsigned time = (run_end - run_begin).usecval();
            if (unlikely(balancer->_verbose))
                click_chatter("Imbalance at run %d : %f-%f %f-%f, square %f, m %f, in %d usec",run,imbalance_o,imbalance_u,overload_allowed,underload_allowed,square_imbalance, m, time);


            auto &v = (*balancer->_stats)[run - 1];
                v.imbalance += square_imbalance;
                v.count ++;
                v.time += time;


            if (run == max_runs || square_imbalance < target_imbalance) break;

            if (run == 1) {
                overload_allowed = imbalance_o;
                underload_allowed = imbalance_u;
                //bottom_square_imbalance = square_imbalance;
                //min_sq = square_imbalance; //min is the sq for the min allowed, not the min seen sq
                phase = 1; //searching top
            } else {

                if (square_imbalance < min_sq) {
                    oa_min = overload_allowed;
                    ua_min = underload_allowed;
                }

                if (phase == 1) {
                    if (square_imbalance <= last_sq) { //Continue finding top
                        overload_allowed = overload_allowed + overload_allowed / 2;
                        underload_allowed = underload_allowed + underload_allowed / 2;
                    } else {
                        phase = 2;
                        top_oa = overload_allowed;
                        top_ua = underload_allowed;
                        m = 0.5;
                        overload_allowed = overload_allowed / 2;
                        underload_allowed = underload_allowed / 2;
                    }
                } else if (phase == 2) { //Searching left inflation
                    if (square_imbalance <= last_sq) {
                        m = (0 + m) / 2; //continue left;
                        if (m < 0.01)
                            break; //Border is the max
                        overload_allowed = bottom_oa + (top_oa - bottom_oa) * m;
                        underload_allowed = bottom_ua + (top_ua - bottom_ua) * m;
                        //Either we still need to descend left
                        //Or we hit the left inflation and will need to descend right afterwards
                    } else if (square_imbalance > last_sq) { //we found a new bottom
                        bottom_oa = overload_allowed;
                        bottom_ua = underload_allowed;
                        bottom_sq = square_imbalance;
                        phase = 3;
                        m = 0.5;
                        overload_allowed = (bottom_oa + top_oa) / 2;
                        underload_allowed = (bottom_ua + top_ua) / 2;
                    }
                } else if (phase == 3) { //Searching right inflation
                    if (square_imbalance <= last_sq) {
                        m = (1 + m) / 2; // continue right
                        if (m > 0.99)
                                break; //Border is the min
                        overload_allowed = bottom_oa + (top_oa - bottom_oa) * m;
                        underload_allowed = bottom_ua + (top_ua - bottom_ua) * m;
                    } else     if (square_imbalance > last_sq) { //we found a new top
                        phase = 2;
                        top_oa = overload_allowed;
                        top_ua = underload_allowed;
                        top_sq = square_imbalance;
                        m = 0.5;
                        overload_allowed = bottom_oa + (top_oa - bottom_oa) * m;
                        underload_allowed = bottom_ua + (top_ua - bottom_ua) * m;
                    }
                }

                if (top_sq == bottom_sq && bottom_sq == square_imbalance && square_imbalance == min_sq) {
                    break;
                }


                /*
                if (square_imbalance <= last_sq) { // Continue in this direction

                } else {
                    //invert_direction, set max according to dir
                    if (dir > 0) {
                        top_oa = overload_allowed;
                        top_ua = underload_allowed;
                        dir = -1;
                        overload_allowed = (bottom_oa + overload_allowed) / 2;
                        underload_allowed = (bottom_ua + underload_allowed) / 2;
                    } else {
                        bottom_oa = overload_allowed;
                        bottom_ua = underload_allowed;
                        dir = 1;
                        overload_allowed = (top_oa + overload_allowed) / 2;
                        underload_allowed = (top_ua + underload_allowed) / 2;
                        bottom_square_imbalance = square_imbalance;
                    }
                    if (bottom_square_imbalance)
                }*/
                /*if (square_imbalance == last_sq)
                    break;*/


                /*if (square_imbalance > last_sq) {
                    dir = -dir;
                    n_change++;
                    if (nchange > 2)
                } else if (square_imbalance < last_sq) {
                    n_change = 0;
                }
                overload_allowed = overload_allowed + dir * (overload_allowed / 2);
                underload_allowed = underload_allowed + dir * (underload_allowed / 2);*/
            }
            if (unlikely(balancer->_verbose > 2))
                click_chatter("Phase %d", phase);

            run++;
            last_sq = square_imbalance;
            if (run == max_runs) {
                if (square_imbalance != min_sq) {
                    overload_allowed = oa_min;
                    underload_allowed = ua_min;
                } else
                    break;
            }

        }


/*
 * Not sure about the gradient
        nlopt::opt opt(nlopt::LD_MMA, 2);
        std::vector<double> lb(2);
        lb[0] = -HUGE_VAL; lb[1] = 0;
        opt.set_lower_bounds(lb);
        opt.set_min_objective(myfunc, NULL);
        my_constraint_data data[2] = { {2,0}, {-1,1} };
        opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
        opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);
        opt.set_xtol_rel(1e-4);
        std::vector<double> x(2);
        x[0] = 1.234; x[1] = 5.678;
        double minf;

        try{
            nlopt::result result = opt.optimize(x, minf);
            std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
                << std::setprecision(10) << minf << std::endl;
        }
        catch(std::exception &e) {
            std::cout << "nlopt failed: " << e.what() << std::endl;
        }*/
        /*
        lsint weights[] = {10, 60, 30, 40, 30, 20, 20, 2};
        lsint values[] = {1, 10, 15, 40, 60, 90, 100, 15};
        lsint knapsackBound = 102;

        // Declares the optimization model.
        LocalSolver localsolver;
        LSModel model = localsolver.getModel();

        // 0-1 decisions
        LSExpression x[8];
        for (int i = 0; i < 8; i++)
            x[i] = model.boolVar();

        // knapsackWeight <- 10*x0 + 60*x1 + 30*x2 + 40*x3 + 30*x4 + 20*x5 + 20*x6 + 2*x7;
        LSExpression knapsackWeight = model.sum();
        for (int i = 0; i < 8; i++)
            knapsackWeight += weights[i]*x[i];

        // knapsackWeight <= knapsackBound;
        model.constraint(knapsackWeight <= knapsackBound);

        // knapsackValue <- 1*x0 + 10*x1 + 15*x2 + 40*x3 + 60*x4 + 90*x5 + 100*x6 + 15*x7;
        LSExpression knapsackValue = model.sum();
        for (int i = 0; i < 8; i++)
            knapsackValue += values[i]*x[i];

        // maximize knapsackValue;
        model.maximize(knapsackValue);

        // close model, then solve
        model.close();

        // Parameterizes the solver.
        localsolver.getParam().setTimeLimit(1);
        localsolver.solve();*/

    /*    glp_prob *lp;
        lp = glp_create_prob();
        glp_set_prob_name(mip, "sample");
        //int x[target.size()]
        int* x = transfer.data();
        glp_set_obj_dir(lp, GLP_MIN);

        //glp_add_rows(lp,buckets_load.size());

        glp_add_cols(lp,buckets_load.size());
        for (int i = 0; i < buckets_load.size(); i++) {
            //glp_set_col_name(mip, 1, "x1");
            // glp_set_col_bnds(mip, 1, GLP_DB, 0.0, 40.0);
             // glp_set_obj_coef(mip, 1, 1.0);
              //glp_set_col_name(lp, 1, "b1");
              glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
              glp_set_obj_coef(lp, i, buckets_load.size());
              glp_set_col_kind(lp, i, GLP_IV);
        }

        glp_iocp parm;
          glp_init_iocp(&parm);
          parm.presolve = GLP_ON;
          int err = glp_intopt(mip, &parm);*/

    }

private:
    float min_cost;
};
/*

class Problem
{ public:
   // Vector<int> oid;
   // Vector<int> uid;
    Vector<int> transfer;
    float min_cost;
    //Vector<float> imbalance;
    //float target;
    //int N;

    Problem() : oid(), uid(), min_cost(FLT_MAX) {

    }

    bool computeSol() {
        float newload[N] = {0.0f};
        for (int c = 0; c < N; c++) { //Sum of imbalance for all cores
            newload[transfer[c]] += imbalance[c];
        }
        float imb = 0;
        for (int c = 0; c < N; c++) {
            imb += newload[c] * newload[c];
        }
        if (imb < min_cost) {
            min_cost = imb;
            return true;
        }
        return false;
    }

    bool solve() {
        tryK(0);
    }
*/
    /**
     * Problem of this method:
     * If all sending cores have quite big buckets, then correction may lead to all
     * cores having too low corrected transfer so no bucket will move. Moving some
     * of them would probably be best.
     *
     * One solution is to consider all buckets that fit the correction, then select a subset
     * that min the invariance of the moving set.
     *
     * But that's again a bad minimization
     *
     * So if no buckets were moved, we do a second pass without the correction
     */
  /*  Vector<float> fixedTransfert() {
        //See above
        float newload[N] = {0.0f};
        for (int c = 0; c < N; c++) { //Sum of imbalance for all cores
            newload[transfer[c]] += imbalance[c];
        }
        Vector<float> fT;
        fT.resize(N,1.0);

        //Eg c0 had  0.15 of imbalance //underloaded, newload -0.03
        //   c1 had -0.10 of imbalance //overloaded, newload 0
        //   c2 has -0.08 of imbalance  //overloaded, newload 0

        // c1 and c2 will give all to c0, but c0 will become overloaded by -0.03.
        //Let's share the final imbalance to -0.03/3, so every of those are at -0.01
                        //
        for (int c = 0; c < N; c++) { //C is the receiving core, ie, for each receiving core C
            if (newload[c] < 0) {

                int nsources = 0;
                float imb = 0; //Imbalance of all senders

                for (int j = 0; j < N; j++) {
                    if (j == c)
                        continue;
                    if (transfer[j] == c) {
                        nsources++;
                        imb += imbalance[j]; // += -0.10 += -0.08  -> -0.18
                    }
                }
                if (nsources == 0)
                    continue;
                float avgimb = (imbalance[c] + imb) / (nsources + 1); //Average imbalance Eg 0.15 - 0.18 / 3 = -0.03 / 3 = -0.01



                //If we just give the factor to all sources, then some would be left with more imbalance than others
                // We want all the sources to have exactly the same amount of final imbalance ->
                //   avgimb
                if (avgimb < -0.001) {
                    for (int j = 0; j < N; j++) {
                        if (j == c)
                            continue;
                        if (transfer[j] == c) {

                            // c1 :: -0.01 / -0.10 -> 0.1 -> 0.9
                            // c2 :: -0.01 / -0.08 -> 0.125 -> 0.875
                            fT[j] = 1 - avgimb / imbalance[j];
                        }
                    }
                }

                //Final example : j gives 0.875*-0.2 == -0.175 --> n = -0.025
                //c reveives -0.175 + 0.15 = -0.025

            }
        }
        return fT;

    }

private:
    bool tryK(int i)
    {
        if (i == oid.size()) {
            return computeSol();
        }
        int best = -1;
        for (int j = 0; j < uid.size(); j++) {
                transfer[oid[i]] = uid[j];
                if (tryK(i + 1)) {
                    best = j;
                }
        }
        if (best > -1) {
            transfer[oid[i]] = uid[best];

            return true;
        }
        return false;
    }
};
*/



/**
 * Problem of moving all given buckets to a set of cores, ensuring load balancing.
 *
 * It is solved by simply pushing load of each bucket one by one to the core
 * with the least load. We start with the most loaded buckets first.
 * Not the most optimal solution, but we do not care as we will rebalance
 * shortly.
 */
class BucketMapProblem
{ public:
    Vector<int> transfer; //existing core id for each buckets
    Vector<float> imbalance; //Imbalance for each existing cores (no holes)
    Vector<float> buckets_load; //Load for each bucket

    BucketMapProblem(int nbuckets, int ncpu) : min_cost(FLT_MAX) {
        transfer.resize(nbuckets);
        buckets_load.resize(nbuckets);
        imbalance.resize(ncpu);
    }

    /*This is WAY too long
    bool tryK(int i)
    {

        if (i == buckets_load.size()) { //All buckets have been mapped
            float newload[imbalance.size()] = {0.0f};
            for (int c = 0; c < buckets_load.size(); c++) { //Sum of imbalance for all cores
                newload[transfer[c]] += buckets_load[c];
            }
            float imb = 0;
            for (int c = 0; c < imbalance.size(); c++) {
                float coreload = imbalance[c] + newload[c];
                imb += (coreload) * (coreload);
            }
            if (imb < min_cost) {
                min_cost = imb;
                return true;
            }
            return false;
        }
        int best = -1;
        for (int j = 0; j < imbalance.size(); j++) {
//                click_chatter("Level %d:%d",i,j);
                transfer[i] = j;
                if (tryK(i + 1)) {
                    best = j;
                }
        }
        if (best > -1) {
            transfer[i] = best;

            return true;
        }
        return false;
    }

    void solve() {
        return tryK(0);
    }
    */

    void solve() {
        typedef struct {
            int id;
            float load;
        } bref;
        auto cmp = [](bref left, bref right) { return left.load > right.load; };
        std::priority_queue<bref, std::vector<bref>, decltype(cmp)> q(cmp);
        for(int i = 0; i < buckets_load.size(); i++) {
            float f = buckets_load[i];

            q.push(bref{.id = i,.load =   f});
        }

        auto cmpc = [](bref left, bref right) { return left.load < right.load; };
        std::priority_queue<bref, std::vector<bref>, decltype(cmp)> cores(cmp);
        for(int i = 0; i < imbalance.size(); i++) {
            float f = imbalance[i];
            //click_chatter("Core %d should receive %f load",i,f);
            cores.push(bref{.id = i,.load = - f});//negative of imbalance, so we should reach a nice 0 everywhere by adding some load
        }


        while (!q.empty()) {
            bref t = q.top();
            q.pop();
            bref c = cores.top();
            cores.pop();
            transfer[t.id] = c.id;
            c.load += t.load;
            cores.push(c);
        }

    }

private:
    float min_cost;
};





#endif
