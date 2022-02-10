# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 16:43:37 2022

@author: drc9
"""

import numpy as np
import os
import gurobipy as gp

folder_name = "test_instance"
alpha_S = 0 # Importance given to student time preference
alpha_F = 1 # Importance given to faculty time preference
gurobi_log_file_name = folder_name + "_log.txt"
output_file_name = folder_name + "_timetable.txt"

input_set_T_D_param_ct_TD = np.genfromtxt(os.path.join(os.getcwd(),folder_name,"1_set_T_D_param_ct_TD.csv"),delimiter=",",dtype=str)
input_set_I = np.genfromtxt(os.path.join(os.getcwd(),folder_name,"2_set_I.csv"),delimiter=",",dtype=str)
input_set_J_Jdash_param_m_pF = np.genfromtxt(os.path.join(os.getcwd(),folder_name,"3_set_J_Jdash_param_m_pF.csv"),delimiter=",",dtype=str)
input_set_K_param_s_ep_pS = np.genfromtxt(os.path.join(os.getcwd(),folder_name,"4_set_K_param_s_ep_pS.csv"),delimiter=",",dtype=str)
input_a = np.genfromtxt(os.path.join(os.getcwd(),folder_name,"5_a.csv"),delimiter=",",dtype=str)
input_cc = np.genfromtxt(os.path.join(os.getcwd(),folder_name,"6_cc.csv"),delimiter=",",dtype=str)
input_cf = np.genfromtxt(os.path.join(os.getcwd(),folder_name,"7_cf.csv"),delimiter=",",dtype=str)


##############################################################################
##############################################################################

def get_set_T_D_param_ct_TD(input_set_T_D_param_ct_TD):
    if len(input_set_T_D_param_ct_TD[0,:]) != 4:
        print("ERROR: 1_set_T_D_param_ct_TD.csv should have exactly 4 columns.")
        return 0,0,0,0
    else:
        set_T = input_set_T_D_param_ct_TD[1:,0]
        set_D = np.unique(input_set_T_D_param_ct_TD[1:,1])
        days = (input_set_T_D_param_ct_TD[1:,1])
        start_times = (input_set_T_D_param_ct_TD[1:,2]).astype(int)
        end_times = (input_set_T_D_param_ct_TD[1:,3]).astype(int)
        ct = np.zeros((len(set_T),len(set_T)))
        TD = 1e+6*np.ones((len(set_T),len(set_D)))
        for t in range(len(set_T)):
            for tdash in range(len(set_T)):
                condition1 = (start_times[tdash] >= start_times[t]) and (start_times[tdash] <= end_times[t])
                condition2 = (end_times[tdash] >= start_times[t]) and (end_times[tdash] <= end_times[t])
                condition3 = (days[tdash] == days[t])
                if ((condition1 or condition2) and condition3):
                    ct[t,tdash] = 1
            arg_day = (np.where(set_D == days[t])[0]).squeeze()
            TD[t,arg_day] = 1  
        return set_T,set_D,ct,TD

def get_set_I(input_set_I):
    if len(input_set_I.shape) != 1:
        print("ERROR: 2_set_I.csv should have exactly 1 column.")
        return 0
    else:
        return (input_set_I[1:])

def get_set_J_Jdash_param_m_pF(input_set_J_Jdash_param_m_pF,set_T):
    if len(input_set_J_Jdash_param_m_pF[0,:]) != 3+len(set_T):
        print("ERROR: for the current version, 3_set_J_Jdash_param_m_pF.csv should have exactly 3+%s columns." %(len(set_T)))
        return 0,0,0,0
    else:
        set_J = input_set_J_Jdash_param_m_pF[1:,0]
        J_dash_index = np.where((input_set_J_Jdash_param_m_pF[1:,1]).astype(int)==1)[0]
        if len(J_dash_index) == 0:
            set_Jdash = np.array([])
        else:
            set_Jdash = set_J[J_dash_index]
        m = (input_set_J_Jdash_param_m_pF[1:,2]).astype(int)
        pF = np.zeros((len(set_J),len(set_T)))
        for t in range(len(set_T)):
            arg_t = (np.where(input_set_J_Jdash_param_m_pF[0,:] == set_T[t])[0]).squeeze()
            pF[:,arg_t-3] = (input_set_J_Jdash_param_m_pF[1:,arg_t]).astype(int)
        return set_J,set_Jdash,m,pF

def get_input_set_K_param_s_ep_pS(input_set_K_param_s_ep_pF,set_I,set_T):
    if len(input_set_K_param_s_ep_pF[0,:]) != (2+len(set_I)+len(set_T)):
        print("ERROR: for the current version, 4_set_J_Jdash_param_m_pF.csv should have exactly 3+%s+%s columns." %(len(set_I),len(set_T)))
        return 0,0,0,0
    else:
        set_K = input_set_K_param_s_ep_pF[1:,0]
        s = (input_set_K_param_s_ep_pF[1:,1]).astype(int)
        ep = np.zeros((len(set_K),len(set_I)))
        for i in range(len(set_I)):
            arg_i = (np.where(input_set_K_param_s_ep_pF[0,:] == set_I[i])[0]).squeeze()
            ep[:,arg_i-2] = (input_set_K_param_s_ep_pF[1:,arg_i]).astype(int)
        pS = np.zeros((len(set_K),len(set_T)))
        for t in range(len(set_T)):
            arg_t = (np.where(input_set_K_param_s_ep_pF[0,:] == set_T[t])[0]).squeeze()
            pS[:,arg_t-2-len(set_I)] = (input_set_K_param_s_ep_pF[1:,arg_t]).astype(int)
        return set_K,s,ep,pS

def get_a(input_a):
    if (len(input_a[0,:]) != (1+len(set_T))) or (len(input_a[:,0]) != (1+len(set_I))):
        print("ERROR: 5_a.csv should have exactly 1+%s rows and 1+%s columns." %(len(set_I),len(set_T)))
        return 0
    else:
        return (input_a[1:,1:]).astype(int)

def get_cc(input_cc,set_I):
    if len(input_cc[0,:]) != 2:
        print("ERROR: 6_cc.csv should have exactly 2 columns.")
        return 0
    else:
        cc = np.zeros((len(set_I),len(set_I)))
        for lines in range(len(input_cc[1:,0])):
            arg_course1 = (np.where(set_I == input_cc[1+lines,0])[0]).squeeze()
            arg_course2 = (np.where(set_I == input_cc[1+lines,1])[0]).squeeze()
            cc[arg_course1,arg_course2] = 1
            cc[arg_course2,arg_course1] = 1
        return cc

def get_cf(input_cf,set_I,set_J):
    if len(input_cf[0,:]) != 2:
        print("ERROR: 7_cf.csv should have exactly 2 columns.")
    else:
        cf = np.zeros((len(set_I),len(set_J)))
        for lines in range(len(input_cf[1:,0])):
            arg_faculty = (np.where(set_J == input_cf[1+lines,0])[0]).squeeze()
            arg_course = (np.where(set_I == input_cf[1+lines,1])[0]).squeeze()
            cf[arg_course,arg_faculty] = 1
        return cf 

def check_input_correctness(set_T,set_D,ct,TD,set_I,set_J,set_Jdash,m,pF,set_K,s,ep,pS,a,cc,cf):
    
    vals = np.zeros(16)
    dictionary = {'a1':set_T, 'a2':set_D, 'a3':ct, 'a4':TD, 'a5':set_I,
                  'a6':set_J, 'a7':set_Jdash, 'a8':m, 'a9':pF, 'a10':set_K,
                  'a11':s, 'a12':ep, 'a13':pS, 'a14':a, 'a15':cc, 'a16':cf}
    q = 0
    for key, value in dictionary.items():
        vals[q] = type(value).__module__ == np.__name__
        q += 1
    
    if np.any(vals):
        return 1
    else:
        return 0

def CSP(alpha_F,alpha_S,gurobi_log_file_name,granular_check):
    
    if granular_check == 0:
        print("ERROR: Atleast one of the input parameters is not defined properly.")
        return 0,0,0
    
    else:
        T = range(len(set_T))
        D = range(len(set_D))
        I = range(len(set_I))
        J = range(len(set_J))
        Jdash = range(len(set_Jdash))
        K = range(len(set_K))
        
        model = gp.Model()
        model.setParam('LogFile',gurobi_log_file_name)
        model.setParam('LogToConsole',0)
        model.setParam('TimeLimit',1200)
        
        x = model.addVars(I,T,vtype=gp.GRB.BINARY,name="x")
        y = model.addVars(J,I,T,vtype=gp.GRB.BINARY,name="y")
        
        obj_1 = alpha_S*gp.quicksum(((pS[k,t]*ep[k,i]/s[k])*x[i,t]) for k in K for i in I for t in T)
        obj_2 = alpha_F*gp.quicksum((pF[j,t]*y[j,i,t]) for j in J for i in I for t in T)
        model.setObjective(obj_1+obj_2, gp.GRB.MAXIMIZE)
        
        model.addConstrs(((gp.quicksum(x[i,t] for t in T) == 1) for i in I), name="eq2")
        model.addConstrs(((gp.quicksum(y[j,i,t] for j in J for t in T) == 1) for i in I), name="eq3")
        model.addConstrs(((x[i,t] + ct[t,tdash]*cc[i,idash]*x[idash,tdash] <= 1) for i in I for idash in I for t in T for tdash in T), name="eq4")
        model.addConstrs(((y[j,i,t] <= x[i,t]) for j in J for i in I for t in T), name="eq5")
        model.addConstrs(((y[j,i,t] <= cf[i,j]) for j in J for i in I for t in T), name="eq6")
        model.addConstrs(((gp.quicksum(y[j,i,t] for i in I) <= 1) for j in J for t in T), name="eq7")
        model.addConstrs(((gp.quicksum(y[j,i,t] for i in I for t in T) <= m[j]) for j in J), name="eq8")
        model.addConstrs(((gp.quicksum(TD[t,d]*y[j,i,t] for i in I for t in T) <= m[j]) for j in Jdash for d in D), name="eq9")
        
        model.optimize()
        
        x_f = np.zeros((len(I),len(T)))
        y_f = np.zeros((len(J),len(I),len(T)))
        for i in I:
            for t in T:
                if x[i,t].X > 0.1:
                    x_f[i,t] = 1
                    for j in J:
                        if y[j,i,t].X > 0.1:
                            y_f[j,i,t] = 1
        
        objval = model.ObjVal
        
        return x_f,y_f,objval
    
def timetable_output(x,y,input_set_T_D_param_ct_TD,output_file_name):
    
    if (type(x).__module__ != np.__name__) or (type(y).__module__ != np.__name__):
        print("ERROR: Atleast one of the input parameters is not defined properly.")
        return 0
    else:
        days_t = (input_set_T_D_param_ct_TD[1:,1])
        start_times = (input_set_T_D_param_ct_TD[1:,2]).astype(int)
        end_times = (input_set_T_D_param_ct_TD[1:,3]).astype(int)
        file = open(output_file_name,"a")
        file.write("\n\n~~~~~~~~ New Timetable Starts ~~~~~~~~\n\n")
        file.write("alpha_S=,%s,alpha_F=,%s\n\n" %(alpha_S,alpha_F))
        for day in range(len(set_D)):
            file.write("Day of week: %s\n" %set_D[day])
            time_index = np.where(days_t == set_D[day])[0]
            start_time_short = start_times[time_index]
            end_time_short = end_times[time_index]
            print(time_index)
            ff = sorted(zip(start_time_short,end_time_short,time_index))
            q = 0
            for a,b,c in ff:
                start_time_short[q] = a
                end_time_short[q] = b
                time_index[q] = c
                q += 1
            for t in time_index:
                file.write("%s-%s," %(start_times[t],end_times[t]))
                for i in np.where(x[:,t]==1)[0]:
                    j = (np.where(y[:,i,t] == 1)[0]).squeeze()
                    file.write("%s:%s," %(set_I[i],set_J[j]))
                file.write("\n")
            file.write("\n")
        file.write("\n~~~~~~~~ Timetable Ends ~~~~~~~~\n\n")
        file.close()
    
        return 0


##############################################################################
##############################################################################

set_T,set_D,ct,TD = get_set_T_D_param_ct_TD(input_set_T_D_param_ct_TD)
set_I = get_set_I(input_set_I)
set_J,set_Jdash,m,pF = get_set_J_Jdash_param_m_pF(input_set_J_Jdash_param_m_pF,set_T)
set_K,s,ep,pS = get_input_set_K_param_s_ep_pS(input_set_K_param_s_ep_pS,set_I,set_T)
a = get_a(input_a)
cc = get_cc(input_cc,set_I)
cf = get_cf(input_cf,set_I,set_J)
granular_check = check_input_correctness(set_T,set_D,ct,TD,set_I,set_J,set_Jdash,m,pF,set_K,s,ep,pS,a,cc,cf)

x_f,y_f,objval = CSP(alpha_F,alpha_S,gurobi_log_file_name,granular_check)
time_table = timetable_output(x_f,y_f,input_set_T_D_param_ct_TD,output_file_name)