
#include "Solver.h"

#include <ctime>
#include <chrono>
#include <iomanip>
//#include <eigen3/Eigen/CholmodSupport>
//#include <gurobi_c++.h>

typedef std::chrono::high_resolution_clock Clock;
//static GRBModel static_model;

void LinearVec::set_label(int label){
    this->label = label;
}

void LinearVec::set_ab(double a, double b){
    this->a = a;
    this->b = b;
}

void LinearVec::clear(){
    ts.clear();
    ws.clear();
}

void LinearVec::push_back(int t, double w){
    ts.push_back(t);
    ws.push_back(w);
}

void QuadraticVec::push_back(int i,int j,double val){
    ts.push_back(triple(i,j,val));
}

void QuadraticVec::clear(){
    ts.clear();
}

void QuadraticVec::set_b(double b){
    this->b = b;
}

void Solution_Struct::init(int n){
    solveval.reserve(n);
    Statue = 0;
}


int optwrapper(vector<double>&para,
               vector<double>&lowerbound,
               vector<double>&upperbound,
               void *funcPara,
               double (*optfunc)(const vector<double>&,vector<double>&,void *),
               double (*constraintfunc)(const vector<double>&,vector<double>&,void *),
               double tor,
               int maxIter,
               double &finEnergy
               ){

    nlopt::result result;
    try{
        int numOfPara = para.size();
        //nlopt::opt myopt(nlopt::LD_CCSAQ,uint(numOfPara));
        //nlopt::opt myopt(nlopt::LD_SLSQP,uint(numOfPara));
        nlopt::opt myopt(nlopt::LD_LBFGS,uint(numOfPara));
        myopt.set_ftol_rel(tor);
        //myopt.set_ftol_abs(tor);

        //myopt.set_xtol_abs(1e-6);
        //myopt.set_xtol_rel(1e-7);

        myopt.set_maxeval(maxIter);

        //myopt.set_initial_step(0.001);
        myopt.set_lower_bounds(lowerbound);
        myopt.set_upper_bounds(upperbound);


        myopt.set_min_objective(optfunc,funcPara);
        //myopt.set_precond_min_objective(optfuncModify, pre, &funcPara);
        //myopt.set_min_objective(optfunc, &funcPara);
        myopt.add_equality_constraint(constraintfunc,funcPara,tor);



        result = myopt.optimize(para, finEnergy);
        std::cout << "Obj: "<< std::setprecision(10) << finEnergy << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }



    return result;


}

int Solver::nloptwrapper(vector<double>&lowerbound,
                         vector<double>&upperbound,
                         nlopt::vfunc optfunc,
                         void *funcPara,
                         vector<NLOptConstraint>&nl_constraints,
                         double tor,
                         int maxIter,
                         Solution_Struct &sol
                         ){

    nlopt::result result;
    sol.Statue=0;

    vector<double>tmp_grad(0);
    sol.init_energy = optfunc(sol.solveval,tmp_grad,funcPara);
    try{
        int numOfPara = lowerbound.size();
        //nlopt::opt myopt(nlopt::LD_CCSAQ,uint(numOfPara));
        nlopt::opt myopt(nlopt::LD_SLSQP,uint(numOfPara));
        //nlopt::opt myopt(nlopt::LD_VAR2,uint(numOfPara));
        myopt.set_ftol_rel(tor);
        //myopt.set_ftol_abs(tor);

        //myopt.set_xtol_abs(1e-6);
        //myopt.set_xtol_rel(1e-7);

        myopt.set_maxeval(maxIter);

        //myopt.set_initial_step(0.001);
        myopt.set_lower_bounds(lowerbound);
        myopt.set_upper_bounds(upperbound);


        myopt.set_min_objective(optfunc,funcPara);
        //myopt.set_precond_min_objective(optfuncModify, pre, &funcPara);
        //myopt.set_min_objective(optfunc, &funcPara);

        for(auto &a:nl_constraints){
            if(a.mt==CM_EQUAL){
                myopt.add_equality_constraint(a.constraintfunc,a.data,tor);
            }else {
                myopt.add_inequality_constraint(a.constraintfunc,a.data,tor);
            }
        }


        auto t1 = Clock::now();

        result = myopt.optimize(sol.solveval, sol.energy);

        auto t2 = Clock::now();

        cout << "nlopt time: " << (sol.time = std::chrono::nanoseconds(t2 - t1).count()/1e9) <<endl;




        sol.Statue = (result >= nlopt::SUCCESS);
        cout<<"Statu: "<<result<<endl;
        std::cout << "Obj: "<< std::setprecision(10) << sol.energy << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }



    return result;


}




int Solver::nloptwrapper(vector<double>&lowerbound,
               vector<double>&upperbound,
               nlopt::vfunc optfunc,
               void *funcPara,
               double tor,
               int maxIter,
               Solution_Struct &sol
               ){


    nlopt::result result;
    sol.Statue=0;
    vector<double>tmp_grad(0);

    sol.init_energy = optfunc(sol.solveval,tmp_grad,funcPara);
    try{
        int numOfPara = lowerbound.size();
        //nlopt::opt myopt(nlopt::LD_VAR1,uint(numOfPara));
        //nlopt::opt myopt(nlopt::LD_CCSAQ,uint(numOfPara));
        //nlopt::opt myopt(nlopt::LD_SLSQP,uint(numOfPara));
        nlopt::opt myopt(nlopt::LD_LBFGS,uint(numOfPara));
        myopt.set_ftol_rel(tor);
        //myopt.set_ftol_abs(tor);

        //myopt.set_xtol_abs(1e-6);
        //myopt.set_xtol_rel(1e-7);
        myopt.set_maxeval(maxIter);

        //myopt.set_initial_step(0.001);
        myopt.set_lower_bounds(lowerbound);
        myopt.set_upper_bounds(upperbound);

        myopt.set_min_objective(optfunc,funcPara);
        //myopt.set_precond_min_objective(optfuncModify, pre, &funcPara);
        //myopt.set_min_objective(optfunc, &funcPara);


        auto t1 = Clock::now();

        result = myopt.optimize(sol.solveval, sol.energy);

        auto t2 = Clock::now();

        cout << "nlopt time: " << (sol.time = std::chrono::nanoseconds(t2 - t1).count()/1e9) <<endl;




        sol.Statue = (result >= nlopt::SUCCESS);
        cout<<"Statu: "<<result<<endl;
        std::cout << "Obj: "<< std::setprecision(10) << sol.init_energy << " -> " <<sol.energy << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }



    return result;


}









//int solveQuadraticProgramming_Core(GRBModel &model, vector<GRBVar>&vars, Solution_Struct &sol, bool suppressinfo = false){

//    int n = vars.size();

//    try{

//        auto t1 = Clock::now();

//        model.optimize();

//        auto t2 = Clock::now();

//        if(!suppressinfo)cout << "Total opt time: " << (sol.time = std::chrono::nanoseconds(t2 - t1).count()/1e9) <<endl<<endl;
//        if(!suppressinfo)cout << "Status: " << model.get(GRB_IntAttr_Status) << endl;

//        //cout<<"0"<<endl;
//        if (sol.Statue = (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)) {
//            sol.solveval.resize(n);
//            //cout<<"1"<<endl;
//            sol.energy = model.get(GRB_DoubleAttr_ObjVal);
//            if(!suppressinfo)cout << "Obj: " << sol.energy << endl;
//            //cout<<"2"<<endl;
//            for(int i=0;i<n;++i)sol.solveval[i] = vars[i].get(GRB_DoubleAttr_X);
//            //for(int i=0;i<n;++i)cout<<sol.solveval[i]<<' ';cout<<endl;

//        }

//    }catch (GRBException e) {
//        cout << "Error code = " << e.getErrorCode() << endl;
//        cout << e.getMessage() << endl;
//    } catch (...) {
//        cout << "Exception during optimization" << endl;
//    }

//    return model.get(GRB_IntAttr_Status);
//    //return time;

//}
//int Solver::solveQuadraticProgramming(vector<triple> &M, vector<triple> &Ab, int n, vector<double>&solveval){


//    vector<GRBVar>vars(n);

//    try{

//        GRBEnv env = GRBEnv();
//        //cout << "optimize 0" << endl;
//        GRBModel model = GRBModel(env);

//        model.getEnv().set(GRB_IntParam_OutputFlag, 0);



//        for(int i=0;i<n;++i)vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,to_string(i));


//        //cout << "optimize 1" << endl;
//        GRBQuadExpr obj = 0;

//        for(int i=0;i<M.size();++i){

//            obj += vars[M[i].i] * vars[M[i].j] * M[i].val;

//        }
//        // cout << "optimize 2 " <<codeer<< endl;
//        model.setObjective(obj, GRB_MINIMIZE);


//        //model.update();

//        for(int i=0;i<Ab.size();++i){
//            model.addConstr(vars[Ab[i].i] - vars[Ab[i].j] >= Ab[i].val);
//        }



//        auto t1 = Clock::now();

//        model.optimize();

//        auto t2 = Clock::now();

//        cout << "Total opt time: " << std::chrono::nanoseconds(t2 - t1).count()/1e9<< endl<< endl;
//        //cout << "Status: " << model.get(GRB_IntAttr_Status) << endl;
//        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
//            solveval.resize(n);
//            cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
//            for(int i=0;i<n;++i)solveval[i] = vars[i].get(GRB_DoubleAttr_X) ;
//        }


//    }catch (GRBException e) {
//        cout << "Error code = " << e.getErrorCode() << endl;
//        cout << e.getMessage() << endl;
//    } catch (...) {
//        cout << "Exception during optimization" << endl;
//    }


//	return 1;
//}


//int Solver::solveQuadraticProgramming(vector<triple> &M, vector<LinearVec> &Ab, int n, vector<double>&solveval){

//    vector<GRBVar>vars(n);

//    try{

//        GRBEnv env = GRBEnv();
//        //cout << "optimize 0" << endl;
//        GRBModel model = GRBModel(env);

//        model.getEnv().set(GRB_IntParam_OutputFlag, 0);

//        for(int i=0;i<n;++i)vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,to_string(i));

//        //cout << "optimize 1" << endl;
//        GRBQuadExpr obj = 0;

//        for(int i=0;i<M.size();++i){

//            obj += vars[M[i].i] * vars[M[i].j] * M[i].val;

//        }
//        // cout << "optimize 2 " <<codeer<< endl;
//        model.setObjective(obj, GRB_MINIMIZE);


//        //model.update();

//        for(int i=0;i<Ab.size();++i){
//            auto &a = Ab[i];
//            GRBLinExpr cons = 0;
//            for(int i=0;i<a.ts.size();++i)cons += vars[a.ts[i]] * a.ws[i];
//            if(a.label == 0)model.addConstr(cons == a.b);
//            else if(a.label == -1)model.addConstr(cons <= a.b);
//            else if(a.label == 1)model.addConstr(cons >= a.b);
//        }



//        auto t1 = Clock::now();

//        model.optimize();

//        auto t2 = Clock::now();

//        cout << "Total opt time: " << std::chrono::nanoseconds(t2 - t1).count()/1e9<< endl<< endl;
//        //cout << "Status: " << model.get(GRB_IntAttr_Status) << endl;
//        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
//            solveval.resize(n);
//            cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
//            for(int i=0;i<n;++i)solveval[i] = vars[i].get(GRB_DoubleAttr_X) ;
//        }

//    }catch (GRBException e) {
//        cout << "Error code = " << e.getErrorCode() << endl;
//        cout << e.getMessage() << endl;
//    } catch (...) {
//        cout << "Exception during optimization" << endl;
//    }
//	return 1;
//}

//int Solver::solveQuadraticProgramming(arma::mat &M, vector<LinearVec> &Ab, int n, Solution_Struct &sol){

//    vector<GRBVar>vars(n);
//    int re;

//    sol.init(n);

//    try{

//        GRBEnv env = GRBEnv();
//        //cout << "optimize 0" << endl;
//        GRBModel model = GRBModel(env);

//        model.getEnv().set(GRB_IntParam_OutputFlag, 0);

//        for(int i=0;i<n;++i)vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,to_string(i));

//        //cout << "optimize 1" << endl;
//        GRBQuadExpr obj = 0;

//        int mm = M.n_rows;
//        for(int i=0;i<mm;++i)for(int j=0;j<mm;++j){

//            obj += vars[i] * vars[j] * M(i,j);

//        }
//        // cout << "optimize 2 " <<codeer<< endl;
//        model.setObjective(obj, GRB_MINIMIZE);


//        int ec = 0, iec = 0;
//        for(int i=0;i<Ab.size();++i){
//            auto &a = Ab[i];
//            GRBLinExpr cons = 0;
//            for(int i=0;i<a.ts.size();++i)cons += vars[a.ts[i]] * a.ws[i];
//            //cout<<a.ts.size()<<endl;
//            if(a.a == a.b){model.addConstr(cons == a.b);ec++;}
//            else {model.addConstr(cons >= a.a);model.addConstr(cons <= a.b);iec+=2;}

//        }
//        cout<<ec<<' '<<iec<<endl;


//        cout<<"Begin Solve"<<endl;
//        re = solveQuadraticProgramming_Core(model, vars, sol);

//        if(0){
//            for(int i=0;i<Ab.size();++i){
//                auto &a = Ab[i];
//                double cons = 0;
//                for(int i=0;i<a.ts.size();++i)cons += vars[a.ts[i]].get(GRB_DoubleAttr_X) * a.ws[i];
//                cout<<cons<<endl;
//            }
//        }



//    }catch (GRBException e) {
//        cout << "Error code = " << e.getErrorCode() << endl;
//        cout << e.getMessage() << endl;
//    } catch (...) {
//        cout << "Exception during optimization" << endl;
//    }

//    return re;


//}





//int Solver::solveLinearLS(vector<triple> &M, vector<double> &b, int n, vector<double>&solveval){

//	return -1;


//}



//int Solver::solveQCQP(arma::mat &M, vector<LinearVec> &Ab, vector<QuadraticVec> &Cb, int n, Solution_Struct &sol){
//    vector<GRBVar>vars(n);
//    int re;

//    sol.init(n);

//    try{

//        GRBEnv env = GRBEnv();
//        //cout << "optimize 0" << endl;
//        GRBModel model = GRBModel(env);

//        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
//        //model.getEnv().set(GRB_INT_PAR_DUALREDUCTIONS, 0);

//        for(int i=0;i<n;++i)vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,to_string(i));

//        //cout << "optimize 1" << endl;
//        GRBQuadExpr obj = 0;

//        int mm = M.n_rows;
//        for(int i=0;i<mm;++i)for(int j=0;j<mm;++j){

//            obj += vars[i] * vars[j] * M(i,j);

//        }
//        // cout << "optimize 2 " <<codeer<< endl;
//        //model.setObjective(obj, GRB_MINIMIZE);
//        model.set(GRB_IntParam_DualReductions, 0);
//        model.setObjective(obj, GRB_MAXIMIZE);
//        //cout<<"model.set(GRB_INT_PAR_DUALREDUCTIONS,0): "<<model.get(GRB_IntA)<<endl;


//        int ec = 0, iec = 0;
//        for(int i=0;i<Ab.size();++i){
//            auto &a = Ab[i];
//            GRBLinExpr cons = 0;
//            for(int i=0;i<a.ts.size();++i)cons += vars[a.ts[i]] * a.ws[i];
//            model.addConstr(cons <= a.a);
//            //cout<<a.ts.size()<<endl;
//            if(a.mt == CM_EQUAL){model.addConstr(cons == a.a);ec++;}
//            else if(a.mt == CM_LESS){model.addConstr(cons <= a.b);iec++;}
//            else if(a.mt == CM_GREATER){model.addConstr(cons >= a.a);iec++;}
//            else if(a.mt == CM_RANGE){model.addConstr(cons >= a.a);model.addConstr(cons <= a.b);iec+=2;}
//        }
//        cout<<"equalc & inequalc: "<<ec<<' '<<iec<<endl;
//        for(int i=0;i<Cb.size();++i){
//            QuadraticVec &a = Cb[i];
//            GRBQuadExpr cons = 0;
//            for(int i=0;i<a.ts.size();++i)cons += vars[a.ts[i].i] * vars[a.ts[i].j] * a.ts[i].val;
//            if(a.mt == CM_EQUAL){model.addQConstr(cons == a.b);ec++;}
//            else if(a.mt == CM_LESS){model.addQConstr(cons <= a.b);iec++;}
//            else if(a.mt == CM_GREATER){model.addQConstr(cons >= a.b);iec++;}
//        }
//        cout<<"equalc & inequalc: "<<ec<<' '<<iec<<endl;



//        cout<<"Begin Solve"<<endl;
//        re = solveQuadraticProgramming_Core(model, vars, sol);

//        if(0){
//            for(int i=0;i<Ab.size();++i){
//                auto &a = Ab[i];
//                double cons = 0;
//                for(int i=0;i<a.ts.size();++i)cons += vars[a.ts[i]].get(GRB_DoubleAttr_X) * a.ws[i];
//                cout<<cons<<endl;
//            }
//        }



//    }catch (GRBException e) {
//        cout << "Error code = " << e.getErrorCode() << endl;
//        cout << e.getMessage() << endl;
//    } catch (...) {
//        cout << "Exception during optimization" << endl;
//    }

//    return re;


//}


//int Solver::solveQP_multiupdata_main( Solution_Struct &sol, vector<LinearVec> &Ab, int n, bool suppressinfo){


//    vector<GRBVar>vars(n);
//    int re;

//    sol.init(n);

//    try{


//        //cout << "optimize 0" << endl;
//        GRBModel model(objmodel);

//        //cout << "optimize 1" << endl;
//        GRBVar *pVar = model.getVars();
//        for(int i=0;i<n;++i){
//            vars[i] = pVar[i];
//        }

//        //model.update();

//        //cout<<"Begin Cons"<<endl;

//        int ec = 0, iec = 0;
//        for(int i=0;i<Ab.size();++i){
//            auto &a = Ab[i];
//            GRBLinExpr cons = 0;
//            for(int i=0;i<a.ts.size();++i)cons += vars[a.ts[i]] * a.ws[i];
//            //cout<<a.ts.size()<<endl;
//            if(a.a == a.b){model.addConstr(cons == a.b);ec++;}
//            else {model.addConstr(cons >= a.a);model.addConstr(cons <= a.b);iec+=2;}

//        }
//        //cout<<ec<<' '<<iec<<endl;


//        //cout<<"Begin Solve"<<endl;
//        re = solveQuadraticProgramming_Core(model, vars, sol, suppressinfo);

//        if(0){
//            for(int i=0;i<Ab.size();++i){
//                auto &a = Ab[i];
//                double cons = 0;
//                for(int i=0;i<a.ts.size();++i)cons += vars[a.ts[i]].get(GRB_DoubleAttr_X) * a.ws[i];
//                cout<<cons<<endl;
//            }
//        }



//    }catch (GRBException e) {
//        cout << "Error code = " << e.getErrorCode() << endl;
//        cout << e.getMessage() << endl;
//    } catch (...) {
//        cout << "Exception during optimization" << endl;
//    }

//    return re;



//}


//int Solver::solveQP_multiupdata_init(arma::mat &M, int n){


//    c_vars.resize(n);
//    for(int i=0;i<n;++i)c_vars[i] = objmodel.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,to_string(i));

//    //cout << "optimize 1" << endl;
//    GRBQuadExpr obj = 0;

//    int mm = M.n_rows;
//    for(int i=0;i<mm;++i)for(int j=0;j<mm;++j){

//        obj += c_vars[i] * c_vars[j] * M(i,j);

//    }
//    // cout << "optimize 2 " <<codeer<< endl;
//    objmodel.setObjective(obj, GRB_MINIMIZE);

//    objmodel.update();
//    cout<<"solveQP_multiupdata_init"<<endl;





//}

