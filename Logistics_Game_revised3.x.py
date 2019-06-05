# -*- coding: utf-8 -*-
from gurobipy import Model, quicksum, GRB
import pandas as pd


DEMAND_FILE = 'demand.csv'

def make_LogiModel(I,J,K,T,hC,hD,b,a,FO,FC,Mk,fik,TCCi,TCDij,djt):
    '''
    Sets
        I          : 流通センター (i)
        J          : 販社 (j)
        T          : 計画期 (t,u)
        K          : センタータイプ (k)
    Parameters
        TCCi       : 工場からセンターiへの輸送費用 (万円/個)
        TCDij      : センターiから販社jへの輸送費用 (万円/個)
        hC         : センター在庫費用 (万円/個)
        hD         : 販社在庫費用 (万円/個)
        b          : 品切れ費用 (万円/個)
        a          : 廃棄費用 (万円/個)
        fik        : タイプkのセンターiの施設費用 (万円/期)
        FO         : センター新規契約費用 (万円/回)
        FC         : センター解約費用 (万円/回)
        Mk         : タイプkのセンターの取り扱い上限量 (個)
        djt        : 販社jのt期の需要量 (個)
    Variables
        X[i,t]     : t期のセンターiへの輸送量
        x[i,j,t]   : t期のセンターiから販社jへの配送量
        y[i,k,t]   : t期にタイプkのセンターiを利用するか否か
        z[i,k,t]   : t期にタイプkのセンターiを新規契約するか否か
        w[i,k,t]   : t期にタイプkのセンターiを解約するか否か
        sC[i,t]    : t期のセンターiの在庫量
        sD[j,t]    : t期の販社jの在庫量
        xi[j,t]    : t期の販社jの品切れ量
        theta[i,t] : t期のセンターiの廃棄量
    '''
    model = Model("Logis")
    
    X,x,y,z,w,sC,sD,xi,theta = {},{},{},{},{},{},{},{},{}

    for t in T:
        for i in I:
            X[i,t]     = model.addVar(vtype=GRB.INTEGER, name="X[{},{}]".format(i,t))
            sC[i,t]    = model.addVar(vtype=GRB.INTEGER, name="sC[{},{}]".format(i,t))
            theta[i,t] = model.addVar(vtype=GRB.INTEGER, name="theta[{},{}]".format(i,t))
            for k in K:
                y[i,k,t] = model.addVar(vtype=GRB.BINARY, name="y[{},{},{}]".format(i,k,t))
                z[i,k,t] = model.addVar(vtype=GRB.BINARY, name="z[{},{},{}]".format(i,k,t))
                w[i,k,t] = model.addVar(vtype=GRB.BINARY, name="w[{},{},{}]".format(i,k,t))
            for j in J:
                x[i,j,t] = model.addVar(vtype=GRB.INTEGER, name="x[{},{},{}]".format(i,j,t))
        for j in J:
            sD[j,t] = model.addVar(vtype=GRB.INTEGER, name="sD[{},{}]".format(j,t))
            xi[j,t] = model.addVar(vtype=GRB.INTEGER, name="xi[{},{}]".format(j,t))
    for i in I:
        sC[i,0] = model.addVar(vtype=GRB.INTEGER, ub=0, name="sC[{},{}]".format(i,0))
        for k in K:
            y[i,k,0] = model.addVar(vtype=GRB.BINARY, ub=0, name="y[{},{},{}]".format(i,k,0))
    for j in J:
        sD[j,0] = model.addVar(vtype=GRB.INTEGER, ub=0, name="sD[{},{}]".format(j,0))

    model.update()

    Only, Capa, Flow = {},{},{}
    for i in I:
        for t in T:
            Only[i,t] = model.addConstr(quicksum(y[i,k,t] for k in K) <= 1,
                                        name="Only[{},{}]".format(i,t))
            Capa[i,t] = model.addConstr(X[i,t]+sC[i,t-1] <= quicksum(Mk[k]*y[i,k,t] for k in K),
                                        name="Capa[{},{}]".format(i,t))
            Flow[i,t] = model.addConstr(X[i,t]+sC[i,t-1] == quicksum(x[i,j,t] for j in J)+sC[i,t]+theta[i,t],
                                        name="Flow[{},{}]".format(i,t))

    flow = {}
    for j in J:
        for t in T:
            flow[j,t] = model.addConstr(quicksum(x[i,j,t] for i in I)+sD[j,t-1] == djt[j,t]+sD[j,t]-xi[j,t],
                                        name="flow[{},{}]".format(j,t))

    strong, connect = {},{}
    for i in I:
        for t in T:
            for j in J:
                strong[i,j,t] = model.addConstr(x[i,j,t] <= quicksum(djt[j,t] for t in T)*quicksum(y[i,k,t] for k in K),
                                                name="Strong[{},{},{}]".format(i,j,t))
            for k in K:
                connect[i,k,t] = model.addConstr(y[i,k,t] - y[i,k,t-1] == z[i,k,t] - w[i,k,t],
                                                 name="Connect[{},{},{}]".format(i,k,t))

    model.setObjective(  quicksum(fik[i,k]*y[i,k,t] for i in I for k in K for t in T)
                       + quicksum(TCCi[i]*X[i,t] for i in I for t in T)
                       + quicksum(TCDij[i,j]*x[i,j,t] for i in I for j in J for t in T)
                       + quicksum(hC*sC[i,t] for i in I for t in T)
                       + quicksum(hD*sD[j,t] for j in J for t in T)
                       + quicksum(b*xi[j,t] for j in J for t in T)
                       + quicksum(FO*z[i,k,t] for i in I for k in K for t in T)
                       + quicksum(FC*w[i,k,t] for i in I for k in K for t in T)
                       + quicksum(a*theta[i,t] for i in I for t in T)
                       ,GRB.MINIMIZE)

    model.__data = X,x,y,z,w,sC,sD,xi,theta
    return model

def print_result(model):
    X,x,y,z,w,sC,sD,xi,theta = model.__data
    TermObj = {t:   quicksum(fik[i,k]*y[i,k,t].X for i in I for k in K)
                  + quicksum(TCCi[i]*X[i,t].X for i in I)
                  + quicksum(TCDij[i,j]*x[i,j,t].X for i in I for j in J)
                  + quicksum(hC*sC[i,t].X for i in I)
                  + quicksum(hD*sD[j,t].X for j in J)
                  + quicksum(b*xi[j,t].X for j in J)
                  + quicksum(FO*z[i,k,t].X for i in I for k in K)
                  + quicksum(FC*w[i,k,t].X for i in I for k in K)
                  + quicksum(a*theta[i,t].X for i in I) for t in T}
    print('-*-*-*-*-*-*- RESULT -*-*-*-*-*-*-')
    print('Total Cost: {0}'.format(model.ObjVal))
    for t in T:
        print('Term {0}'.format(t))
        print('Obj  {0}: {1}'.format(t,TermObj[t]))
        print(' Center [place: type]')
        print('  New   : {0}'.format({int(float(i)):int(float(k)) for i in I for k in K if z[i,k,t].X >= 1}))
        print('  Use   : {0}'.format({int(float(i)):int(float(k)) for i in I for k in K if y[i,k,t].X >= 1}))
        print('  Cancel: {0}'.format({int(float(i)):int(float(k)) for i in I for k in K if w[i,k,t].X >= 1}))
        print(' Transport')
        for i in I:
            if X[i,t].X > 1e-5:
                print('  to Center {0:>2d}: {1:>5s}'.format(int(float(i)),str(X[i,t].X)))
        for j in J:
            for i in I:
                if x[i,j,t].X > 1e-5:
                    print('  to Distributor {0:>2d} from Center {1:>2d}: {2:>5s}'.format(int(float(j)),int(float(i)),str(x[i,j,t].X)))
        print(' Inventories')
        for i in I:
            cStock = False
            if sC[i,t].X > 1e-05:
                cStock = True
                print('  Center {0:>2d}: {1:>5s}'.format(int(float(i)), str(sC[i,t].X)))
        if not cStock:
            print('  Center     : None')
        for j in J:
            dStock = False
            if sD[j,t].X > 1e-05:
                dStock = True
                print('  Distributor {0:>2d}: {1:>5s}'.format(int(float(j)), str(sD[j,t].X)))
        if not dStock:
            print('  Distributor: None')
        print(' Lost Sales')
        for j in J:
            lSales = False
            if xi[j,t].X > 1e-05:
                lSales = True
                print('  Distributor {0:>2d}: {1:>5s}'.format(int(float(j)), str(xi[j,t].X)))
        if not lSales:
            print('  None')
        print(' Disposal')
        for i in I:
            disposal = False
            if theta[i,t].X > 1e-05:
                disposal = True
                print('  Center {0:>2d}: {1:>5s}'.format(int(float(i)), str(theta[i,t].X)))
        if not disposal:
            print('  None')

def read_data():
    hC,hD,b,a,FO,FC = 5,30,200,20,1500,1000

    K = ['1', '2', '3', '4', '5']
    I = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
    Mk = {'5': 30, '2': 100, '1': 120, '4': 40, '3': 80}
    fik = {('6', '1'): 3500, ('2', '5'): 2000, ('7', '3'): 2100, ('8', '3'): 1920, ('6', '3'): 2100, ('8', '2'): 2560, ('7', '1'): 3500, ('1', '2'): 4800, ('11', '4'): 1350, ('11', '1'): 3000, ('5', '4'): 1935, ('10', '1'): 3000, ('4', '2'): 3600, ('10', '3'): 1800, ('8', '5'): 1280, ('10', '5'): 1200, ('5', '2'): 3440, ('9', '1'): 3000, ('11', '2'): 2400, ('8', '4'): 1440, ('3', '2'): 3840, ('6', '5'): 1400, ('1', '1'): 6000, ('7', '2'): 2800, ('9', '3'): 1800, ('9', '5'): 1200, ('7', '4'): 1575, ('6', '4'): 1575, ('10', '2'): 2400, ('7', '5'): 1400, ('2', '4'): 2250, ('3', '1'): 4800, ('2', '1'): 5000, ('3', '4'): 2160, ('9', '2'): 2400, ('5', '1'): 4300, ('4', '3'): 2700, ('2', '3'): 3000, ('4', '4'): 2025, ('9', '4'): 1350, ('5', '3'): 2580, ('4', '1'): 4500, ('2', '2'): 4000, ('3', '3'): 2880, ('11', '5'): 1200, ('5', '5'): 1720, ('8', '1'): 3200, ('3', '5'): 1920, ('1', '4'): 2700, ('4', '5'): 1800, ('6', '2'): 2800, ('10', '4'): 1350, ('1', '5'): 2400, ('11', '3'): 1800, ('1', '3'): 3600}

    J = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
    TCCi = {'5': 13.4, '8': 16.4, '1': 15.2, '4': 12.6, '10': 14.6, '11': 15.4, '7': 16.2, '3': 15.8, '6': 15.8, '2': 14.8, '9': 16.4}
    TCDij = {('6', '4'): 30, ('3', '6'): 28, ('7', '11'): 42, ('1', '7'): 22, ('11', '1'): 18, ('2', '4'): 18, ('11', '10'): 26, ('8', '5'): 34, ('9', '1'): 24, ('9', '3'): 12, ('1', '6'): 18, ('9', '5'): 34, ('2', '5'): 22, ('4', '11'): 16, ('3', '1'): 16, ('2', '1'): 12, ('9', '2'): 14, ('5', '1'): 8, ('3', '9'): 12, ('2', '3'): 8, ('9', '4'): 24, ('5', '3'): 22, ('2', '2'): 0, ('10', '6'): 32, ('5', '7'): 30, ('11', '7'): 42, ('11', '3'): 28, ('7', '8'): 12, ('7', '3'): 22, ('8', '3'): 16, ('8', '2'): 12, ('7', '1'): 22, ('3', '4'): 12, ('10', '7'): 28, ('5', '4'): 16, ('11', '6'): 34, ('1', '8'): 24, ('4', '9'): 24, ('5', '8'): 34, ('1', '9'): 24, ('11', '2'): 30, ('1', '1'): 0, ('7', '2'): 14, ('7', '9'): 20, ('4', '6'): 30, ('2', '11'): 30, ('2', '10'): 14, ('4', '3'): 12, ('4', '8'): 28, ('9', '11'): 40, ('9', '10'): 16, ('4', '1'): 12, ('11', '5'): 12, ('8', '1'): 24, ('10', '11'): 26, ('10', '4'): 10, ('3', '7'): 22, ('7', '10'): 28, ('4', '5'): 16, ('6', '10'): 32, ('8', '8'): 0, ('6', '7'): 16, ('6', '1'): 18, ('11', '8'): 42, ('6', '9'): 32, ('3', '11'): 28, ('6', '3'): 28, ('1', '2'): 12, ('5', '11'): 12, ('11', '4'): 16, ('6', '5'): 22, ('1', '3'): 16, ('10', '5'): 24, ('2', '9'): 14, ('9', '7'): 20, ('2', '8'): 12, ('3', '2'): 8, ('9', '9'): 0, ('10', '10'): 0, ('11', '9'): 40, ('2', '6'): 20, ('11', '11'): 0, ('9', '6'): 32, ('9', '8'): 8, ('5', '9'): 34, ('7', '5'): 30, ('8', '6'): 26, ('2', '7'): 14, ('7', '7'): 0, ('3', '5'): 22, ('5', '6'): 22, ('1', '5'): 8, ('5', '5'): 0, ('10', '3'): 6, ('7', '4'): 30, ('10', '2'): 14, ('8', '7'): 12, ('4', '2'): 18, ('3', '8'): 16, ('1', '10'): 18, ('3', '10'): 6, ('7', '6'): 16, ('5', '2'): 22, ('8', '11'): 42, ('1', '4'): 12, ('6', '11'): 34, ('5', '10'): 24, ('8', '9'): 8, ('4', '7'): 30, ('4', '10'): 10, ('4', '4'): 0, ('3', '3'): 0, ('6', '6'): 0, ('10', '1'): 18, ('8', '10'): 22, ('8', '4'): 28, ('1', '11'): 18, ('10', '9'): 16, ('6', '8'): 26, ('6', '2'): 20, ('10', '8'): 22}

    FILE = pd.read_csv(DEMAND_FILE, header=None)
    NUMBER_OF_TERM = len(FILE.columns)
    T = list(range(1, NUMBER_OF_TERM + 1))
    djt = {(str(i+1),j+1): FILE[j][i] for i in range(11) for j in range(NUMBER_OF_TERM)}
        
    return I,J,K,T,hC,hD,b,a,FO,FC,Mk,fik,TCCi,TCDij,djt

if __name__ == "__main__":
    I,J,K,T,hC,hD,b,a,FO,FC,Mk,fik,TCCi,TCDij,djt = read_data()

    model = make_LogiModel(I,J,K,T,hC,hD,b,a,FO,FC,Mk,fik,TCCi,TCDij,djt)
    model.Params.OutPutFlag = False
    model.optimize()

    print("OptVal:",model.ObjVal)
    print_result(model)
