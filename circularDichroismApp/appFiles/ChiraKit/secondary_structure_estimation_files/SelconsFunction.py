# -*- coding: utf-8 -*-
"""

“The Selcon3 method was originally developed by Sreerama et al (doi: 10.1006/abio.2000.4879) from Colorado State University ”

“The code for this function was originally made in Matlab in 2005 by the research group of
Prof. B.A.Wallace, Birckbeck Collage, London

The MatLab code was updated in 2006-2007 to allow plotting and calculation of
the mean refitted spectrum to the query protein by S.V. Hoffmann, Aarhus University, Denmark

The code was adapted to Python by S.V. Hoffmann, Aarhus University, Denmark in 2021”

This code was later edited by Osvaldo Burastero to allow integration in the ChiraKit online tool, European Molecular Biology Laboratory, Hamburg, 2023

"""
import numpy  as np
import pandas as pd

def SelconsPy(A,F,q, SStruct_labels = ['Alpha-r','Alpha-d', 'Beta-r', 'Beta-d', 'Turns', 'Unord'],Lbl=None):    

    '''
    Input:

        the matrix of reference spectra                 'A' of dimensions m x n (number of wavelengths and proteins)
        the matrix of secondary structure elements      'F' of dimensions l x n (number of secondary structure elements and proteins)
        the query spectrum                              'q' of length m
        the labels for the secondary structure elements 'SStruct_labels' of length l
    '''

    #Selcon2 solutions depend on the RMSD between the protein spectrum and the reconstructed spectrum
    #This is traditionally 0.25 delta epsilon units for CD
    refit_rule=0.25    
    includeProteinsUsedForCalc = 0 #= 1: will add the Protein list to the return string AND if =2 also print to console
    #returns HJ, selcon, selcon2, and selcon3 solutions 
    
    query_spectrum = q

    #******************HJ******************
    #carry out svd on the data
    u, s, vh = np.linalg.svd(A, full_matrices=True)
    #   print ('u-shape: ',u.shape,'s-shape: ', s.shape,'vh-shape: ', vh.shape)
    # This gives u, s and vH such that A = u*S*vh
    # In this case U is 66x66 and vh is 71x71
    # However, s is a verctor with the number of data points (66) singular values. 
    # This has to be made into a matrix, in this case a number of data points x number of ref spectra (66x71)
    S = np.diag(s)
    S = np.zeros((u.shape[0],vh.shape[0])) #66x71 (IR: 200x50)
    s_size = min(u.shape[0],vh.shape[0])
    S[:s_size,:s_size] = np.diag(s) #S[:66,:66] = ... #S[:u.shape[0],:u.shape[0]] = np.diag(s) #S[:66,:66] = ...
    #       print('S-shape: ', S.shape)
    
    # For A = u*S*vh then U*S the basis spectra of the reference dataset 
    # ...and vh contains the least squares coefficients that refit the basis CD spectra to the experimental CD spectra

    #******************HJ5******************
    #Do the matrix multiplication, using the first 5 eigen vectors
    X = F @ np.transpose(vh)[:,:5] @ np.linalg.pinv(S[:5,:5]) @ np.transpose(u[:,:5])
    hj5 = np.matmul(X,q)

    #print('hj5 result')
    SStruct = SStruct_labels
    SSDev   = [x + ' Dev' for x in SStruct_labels]

    format_list = [SStruct[hj5.tolist().index(item)]+'\t'+'{:.4f}' for item in hj5]
    strFormat = '\t'.join(format_list)
    ##print ('hj5 result\t\t  ', strFormat.format(*hj5), 'sum = ', '{:.4f}'.format(hj5.sum()))

    #******************Arrange A (and F) according to similarity to q******************
    #%This is a 1x71 (1x#ref proteins) matrix
    rmsds=np.zeros(A.shape[1]) #rmsds=zeros(1,size(A,2));
    #%calculate index of rmsd values of query to data values
    for i in range(A.shape[1]):  #for i= 1:size(A,2); 
        rmsds[i] = np.sqrt(((q - A[:,i]) ** 2).mean())  #rmsds(i)=rmsd(q,A(:,i)); %SVH 20-9-2006
    
    # Get the indexes of the rmsds sorted so that the smallest rmsds comes first
    ix = np.argsort(rmsds) #[sorted, ix] = sort(rmsds);
    Asort = np.zeros((A.shape[0],A.shape[1])) #Asort = zeros(size(A,1),size(A,2));
    Fsort = np.zeros((F.shape[0],F.shape[1])) #Fsort = zeros(size(F,1),size(F,2));
    #ProtsSort=prots;%  SVH 25/9-06
    
    # Arrange A and F according the the rmsd 
    for i in range(A.shape[1]): #for i= 1:size(A,2);
        Asort[:,i]=A[:,ix[i]] #Asort(:,i)=A(:,ix(i)); 
        Fsort[:,i]=F[:,ix[i]]  #Fsort(:,i)=F(:,ix(i)); 
        #ProtsSort(i)=prots(ix(i));
    #end
    
    #%concatanate q to A
    A=np.c_[q, Asort] #A=[q Asort;];
    #%make an intial guess about ss
    F=np.c_[Fsort[:,0], Fsort] #F=[Fsort(:,1) Fsort;];
    
    #pred.Asort=A;
    #pred.guess = Fsort(:,1);
    
    #******************SELCON1******************
    #%now do selcon1 
    
    #%for %min number of proteins to max  number of proteins
    #%for min number of basis spectra to max number of basis spectra
    #%calculate solution and store it in a matrix to be filtered by sum and
    #%fract rule
    
    selfcon  = 0
    attempts = 0
    while selfcon==0:
        
        #%empty all of the solutions from the variable selection
        solutionsList=[]
        
        #%for various number of cd spectra
        for i in range(2, A.shape[1]): #this is 0 to (71+1)-1    #for i=3:size(A,2)
            #%make the truncated matrix B
            u, s, vh = np.linalg.svd(A[:,:i+1], full_matrices=True)  #[U,S,V] = svd(A(:,1:i),0); #print ('u-shape: ',u.shape,'s-shape: ', s.shape,'vh-shape: ', vh.shape)
            S = np.zeros((u.shape[0],vh.shape[0])) #66x71 #print('S-shape: ', S.shape)
            S[:s.shape[0],:s.shape[0]] = np.diag(s) #S[:66,:66] = ... Note :66 is 0 to 65
            #%calculate the maximum number of basis spectra to use
            maxbasis = 7
            if (i+1 < 7):
                maxbasis = i+1
            #%for various numbers of basis spectra
            for j in range(maxbasis): #this goes from 0 to maxbasis-1 #for j=1:maxbasis
                newSolution = F[:,:i+1] @ np.transpose(vh)[:,:j+1] @ np.linalg.pinv(S[:j+1,:j+1]) @ np.transpose(u[:,:j+1]) @ q #((F(:,1:i)*V(:,1:j) * pinv(S(1:j,1:j)) * U(:,1:j)') *q)
                solutionsList.append(newSolution) #solutions =[solutions ((F(:,1:i)*V(:,1:j) * pinv(S(1:j,1:j)) * U(:,1:j)') *q);];
            #end
        #end
        
        #convert to a matrix (#print('Solutions list')     #print(solutionsList))
        solutions = np.transpose(np.array(solutionsList))
        
        #%now apply the sum and fract rules to the solutions
        #%relaxing as we go to make sure we get a minimum of 1 solution
        
        #%keep track of how far weve iterated along all solutions
        solution_number = 1
        
        #%track all valid solutions and parameters not useful for selcon but might
        #%give better results so try in the future
        valid_solutions = []
        #valid_params   = []
        
        #%store the parameters that have obtained the best solutions for each np
        np_solutions = []
        np_params    = []
        #%init the sum and fract rule that may need to be relaxed later
        sum_rule   = 0.05
        fract_rule = -0.025
        
        #%bool to see if weve found a minimum of one solution after filtering with
        #%the sum and fract rules
        found_one = 0
        while found_one == 0:
            solution_number = 1;    
            #%for each value of Np
            for i in range (2, A.shape[1]): #for i=3:size(A,2)
                #% a way of checking if we have the minimum sum error
                sum_best = 10000
                #% store the best parameter of the basis spectra
                best_basis=0
                #% store the best parameters np and basis for the best solutions
                #best_param=[]
                #% store the best solution for this particular number of np
                best_solution=[]
                #% determine the maximum number of basis spectra for this np
                maxbasis = 7
                if (i+1 < 7):
                    maxbasis=i+1    
                #%for each number of basis spectra
                for j in range(maxbasis): #for j=1:maxbasis
                    #%if the solution satisfies the sum and fraction rule
                    testSum = solutions[:,solution_number-1].sum()
                    testMinFrac = solutions[:,solution_number-1].min()   #print('i,j', i,j,'testSum: ',testSum,' testMinFrac ',testMinFrac)
                    if ((abs(1-testSum) <= sum_rule) and (testMinFrac >=fract_rule)): #if (((abs(1-(sum(solutions(:,solution_number)))))<=sum_rule) && min(solutions(:,solution_number))>=fract_rule);
                        #print('testSum', testSum, 'abs(1-testSum)', abs(1-testSum),'testMinFrac', testMinFrac)
                        #%store the solutions and valid parameters
                        valid_solutions.append(solutions[:,solution_number-1])  ##valid_solutions=[valid_solutions solutions(:,solution_number);];
                        #valid_params.append(np.array([i, j])) #!!!!not sure this is used later!!!!  : valid_params.append([i;j;]) ##valid_params = [valid_params [i;j;];];    
                       
                        #%now keep track of the best solution and its basis number
                        #%for this particular np
                        if abs(1-testSum)<= sum_best: #if (abs(1-(sum(solutions(:,solution_number)))) < sum_best)
                            best_basis = j #best_basis =j;
                            best_solution = solutions[:,solution_number-1] #best_solution = solutions(:,solution_number);
                            sum_best = abs(1-testSum) #sum_best = abs(1-(sum(solutions(:,solution_number))));
                        #end
                    #end
                    #%increment the solution number
                    solution_number = solution_number +1;

                #end %for nj (num of basis spectra, SVH 22/9-2006)
                #%now go through and only include the closest solution to 1 from each value of Np. This is the best_solution
                #% if we have a valid solution for this np then record basis and solution
                if best_basis > 0: #if (best_basis >0)
                    np_params.append(np.array([i, best_basis])) #    np_params= [np_params [i;best_basis;];]; 
                    np_solutions.append(best_solution) #    np_solutions= [np_solutions best_solution;];
                #end
            
            #end%for np (#ref prot used)
            
            if len(valid_solutions) > 0: #if(size(valid_solutions,2)  > 0)
                found_one = 1        #print('used fract_rule', fract_rule, 'sum_rule= ', sum_rule)
            #end
            #%relax
            fract_rule= fract_rule - 0.005
            sum_rule  = sum_rule   + 0.01
        #end %while(found_one==0), SVH 22/9-2006

        #pred.np = np_solutions;
        #pred.np_params = np_params;
        
        #np_solutions is a list of solutions[:,solution_number-1] arrays so create this as a matrix
        np_solutionsMat = np.transpose(np.array(np_solutions))      #print('np_solutionsMat shape: ', np_solutionsMat.shape)
        sel1 = np.mean(np_solutionsMat,axis = 1) #pred.sel1 = mean(np_solutions,2);
        #pred.jon1 = mean(valid_solutions,2);
        
        #%pred.valid = valid_solutions;
        #%pred.valid_params = valid_params;
        
        #disp(['For selfconsist attempt ' num2str(attempts) ' rmsd is ' num2str(rmsd(F(:,1),mean(np_solutions,2)))]) %SVH 22-9-2006
        format_list = [SStruct[sel1.tolist().index(item)]+'\t'+'{:.4f}' for item in sel1]
        strFormat = '\t'.join(format_list)
        ##print ('For selfconsist attempt ', attempts, strFormat.format(*sel1), 'sum = ', '{:.4f}'.format(sel1.sum()))
        
        #%check to see if selfconsistent
        if np.sqrt(((F[:,0] - sel1) ** 2).mean()) <0.0025: #        if((rmsd(F(:,1),mean(np_solutions,2)))<0.0025) %If the new solution is close to the previous, we have selfconsistency, SVH 7/12-17
            selfcon=1;
        #end;
        
        #%can change this bit just to check if using all valid solutions is better
        F[:,0] = sel1 #mean(np_solutions,2); %Make a new guess for the solution, and run through the solution space again, SVH 7/12-17
        attempts =attempts +1;
        if (attempts >=50):
            print('bail*********************************************************')
            selfcon=1
        #end
        
        #selfcon = 1 #dummy end while loop
    #end  %while(selfcon==0)   %while not selfcon

    #******************SELCON2******************
    #%now we can do selcon 2
    #%filer out solutions using the refi errors of the successfuk solutuins
    sel2_solutions = []
    np_param_Sel2  = []  #; %Keep track on valid #Np and #Basis for valid Sel2 solutions, SVH 25/9-06
    #SVH: Note that np_params is a list of 1x2 arrays each array is np and basis
    refits = np.zeros(len(np_params))  #;refits = zeros(1,size(np_params,2));
    
    for i in range(len(np_params)): #for i=1:size(np_params,2)
        np_=np_params[i][0] #np=np_params(1,i);
        basis=np_params[i][1] #basis=np_params(2,i);
        #Note: np is in the range 0 to 71-1, and basis is in the range 0 to 7-1.
        #But when slicing e.g. :3 means 0 to 2, so to include the i and j use +1
        #So use them as indices +1 (old wrong comment: as they are (not +1))
        u, s, vh = np.linalg.svd(A[:,:np_+1], full_matrices=True) #[U,S,V] = svd(A(:,1:np),0);         #Vt=V';
        S = np.zeros((u.shape[0],vh.shape[0])) 
        S[:s.shape[0],:s.shape[0]] = np.diag(s)
        refit = u[:,:basis+1] @ S[:basis+1,:basis+1] @ vh[:basis+1,:] ##refit = U(:,1:basis)*S(1:basis,1:basis)*Vt(1:basis,:);        
        refits[i] = np.sqrt(((refit[:,0] - q) ** 2).mean()) #refits(i)=rmsd(refit(:,1),q); %SVH 22/9-2006        
    #end  

    got_refit=0
    #refit_rule is defined in top of the function
    while got_refit==0:
        for i in range(len(np_params)): #for i=1:size(np_params,2) %size(np_params,2) is the number of selcon1 solutions, SVH 7/12-2017
            if (refits[i] <= refit_rule):
                sel2_solutions.append(np_solutions[i]) #sel2_solutions=[sel2_solutions np_solutions(:,i);];
                np_param_Sel2.append(np_params[i]) #np_param_Sel2=[np_param_Sel2 np_params(:,i);]; %SVH 25/9-06
            #end
        #end; 
        if (len(sel2_solutions) > 0): #if (size(sel2_solutions,2) > 0)
            got_refit=1
        else:
            refit_rule = refit_rule + 0.01 #; %Relax refit rule
        #end
    #end %while(got_refit==0)
    sel2_solutionsMat = np.transpose(np.array(sel2_solutions))
    np_param_Sel2Mat = np.transpose(np.array(np_param_Sel2))
    
    sel2 = np.mean(sel2_solutionsMat,axis = 1)  #pred.sel2 = mean(sel2_solutions,2);
    
    format_list = [SStruct[sel2.tolist().index(item)]+'\t'+'{:.4f}' for item in sel2]
    strFormat = '\t'.join(format_list)
    ##print ('Sel2 solution:\t\t  ', strFormat.format(*sel2), 'sum = ', '{:.4f}'.format(sel2.sum()))
    ##print ('Sel2 solutions finally used RMSD rule:', refit_rule)
    
    #pred.np_param_Sel2=np_param_Sel2;
    #pred.sel2Sol = sel2_solutions; %SVH 7/12-2017
    
    #******************SELCON3******************
    #%and finally we carry out application of the helix rule (selcon3)
    sel3_solutions=[] #; 
    np_param_Sel3=[] #; %SVH 25/9-06
    
    #%************** change these two lines****
    #%not disordered helix
    hjh = hj5[0] #hjh = hj5(1);
    hel = sel2_solutionsMat[0,:] #hel = sel2_solutions(1,:);
    
    #%these two lines need to be uncommented if ad (Aplha Disorted) and (Alpha Regular) assignment is used
    hjh = hj5[0] + hj5[1] #hjh = pred.hj5(1) + pred.hj5(2);
    hel = sel2_solutionsMat[1,:] + sel2_solutionsMat[0,:] #hel = sel2_solutions(2,:) + sel2_solutions(1,:);
    
    #%**************************************************************************
    
    hel_max = np.max(hel)  #hel_max = max(hel);
    hel_min = np.min(hel) #hel_min = min(hel);
    hel_ave = np.mean(hel) #hel_ave = mean(hel);
    
    if (hjh >0.65):
        for i in range(len(sel2_solutions)): #for i=1:size(sel2_solutions,2)
            helix = hel[i]
            if (helix > 0.65):
                sel3_solutions.append(sel2_solutions[i]) #sel3_solutions= [sel3_solutions sel2_solutions(:,i);];
                np_param_Sel3.append(np_param_Sel2[i]) #np_param_Sel3=[np_param_Sel3 np_param_Sel2(:,i);]; %SVH 25/9-06
            #end
        #end
    #end
    
    
    if (hjh <= 0.65 and hjh >= 0.25):
        for i in range(len(sel2_solutions)): #for i=1:size(sel2_solutions,2)
            helix = hel[i]
            if (helix <= (((hjh + hel_max)/2)+0.03) and helix >= (((hjh + hel_max)/2)-0.03)):
                sel3_solutions.append(sel2_solutions[i]) #sel3_solutions= [sel3_solutions sel2_solutions(:,i);];
                np_param_Sel3.append(np_param_Sel2[i]) #np_param_Sel3=[np_param_Sel3 np_param_Sel2(:,i);]; %SVH 25/9-06
            #end
        #end
    #end
    
    if (hjh < 0.25 and hjh >= 0.15):
        for i in range(len(sel2_solutions)): #for i=1:size(sel2_solutions,2)
            helix = hel[i]
            if (helix <= (((hjh + hel_ave)/2)+0.03) and helix >= (((hjh + hel_ave)/2)-0.03)):
                sel3_solutions.append(sel2_solutions[i]) #sel3_solutions= [sel3_solutions sel2_solutions(:,i);];
                np_param_Sel3.append(np_param_Sel2[i]) #np_param_Sel3=[np_param_Sel3 np_param_Sel2(:,i);]; %SVH 25/9-06
            #end
        #end
    #end
    
    if (hjh < 0.15):
        for i in range(len(sel2_solutions)): #for i=1:size(sel2_solutions,2)
            helix = hel[i]
            if (helix <= (((hjh + hel_min)/2)+0.03) and helix >= (((hjh + hel_min)/2)-0.03)):
                sel3_solutions.append(sel2_solutions[i]) #sel3_solutions= [sel3_solutions sel2_solutions(:,i);];
                np_param_Sel3.append(np_param_Sel2[i]) #np_param_Sel3=[np_param_Sel3 np_param_Sel2(:,i);]; %SVH 25/9-06
            #end
        #end
    #end
    
    refit_np_param = [] #; %SVH 25/9-06

    return_List = []
    
    results_df_selcon = None

    if (len(sel3_solutions) > 0): #if(size(sel3_solutions,2) > 0)
        sel3_solutionsMat = np.transpose(np.array(sel3_solutions))
        sel3 = np.mean(sel3_solutionsMat,axis = 1) #  pred.sel3 = mean(sel3_solutions,2);
        
        refit_np_param = np_param_Sel3 #; %SVH 25/9-06
        stDev = np.std(sel3_solutionsMat, axis=1) #std(sel3_solutions,1,2) #;%SVH 8/12/2017 flag=1: divide by n (flag=0 divide by n-1). 2= along second dimension
        #print(stDev)
        #pred.sel3rms=stDev;%SVH 8/12/2017
        #pred.sel3sol=sel3_solutions;
        writeStr = 'There are ' + str(len(sel3_solutions)) + ' Selcon3 solutions'
        return_List.append(writeStr)
        ##print(writeStr)
        #print('There are ', len(sel3_solutions), ' Selcon3 solutions') #disp(['There are ' num2str(size(sel3_solutions,2)) ' Selcon3 solutions']);
        writeStr = '\t Mean     Std'
        return_List.append(writeStr)
        ##print(writeStr)
        #print('\t Mean     Std') #disp ([' Mean'  '     Std']);
        formatSpec = '{:4.1f}'

        for i in range(len(sel3)): #disp ([num2str(100*mean(sel3_solutions,2),formatSpec) percent plusMinus6 num2str(100*stDev,formatSpec) percent;]);
            strFormat = SStruct[i]+'\t'+ formatSpec+' +/-' + formatSpec
            writeStr  = strFormat.format(100*sel3[i], 100*stDev[i])
            return_List.append(writeStr)
            ##print(writeStr)

        SStruct_for_df = SStruct + ['Total']
        sel3_for_df    = np.append(sel3, np.sum(sel3))

        results_df_selcon = pd.DataFrame({'Component_Selcon3':SStruct_for_df,'Percentage':np.round(sel3_for_df*100,1)})

        #print(strFormat.format(100*sel3[i], 100*stDev[i])) 
        writeStr = 'sum =\t' + formatSpec.format(100*np.sum(sel3)) + '%'
        return_List.append(writeStr)
        ##print(writeStr)
        #print('sum =\t', formatSpec.format(100*np.sum(sel3)), '%')   #disp(['sum = ' num2str(sum(100*mean(sel3_solutions,2)),formatSpec) '%']);
        
    else:
        #sel2_solutionsMat and np_param_Sel2Mat already exists, as well as sel2
        sel3 = np.mean(sel2_solutionsMat,axis = 1) # pred.sel3= mean(sel2_solutions,2);
        
        refit_np_param=np_param_Sel2 #; %SVH 25/9-06
        stDev = np.std(sel2_solutionsMat, axis=1) #stDev=std(sel2_solutions,1,2);%SVH 8/12/2017
        #pred.sel2rms=stDev;%SVH 8/12/2017
        #pred.sel2sol=sel2_solutions;
        #fprintf ( 2, 'No Selcon3 solutions!\n' )#;

        #print('No Selcon3 solutions!\n', 'Using the Selcon2 solutions as result') #disp('Using the Selcon2 solutions as result');
        
        writeStr = 'There are ' + str(len(sel2_solutions)) + ' Selcon2 solutions'
        return_List.append(writeStr)
        #print(writeStr)
        #print('There are ', len(sel2_solutions), ' Selcon2 solutions') #disp(['There are ' num2str(size(sel2_solutions,2)) ' Selcon2 solutions']);
        
        writeStr = '\t Mean     Std'
        return_List.append(writeStr)
        #print(writeStr)
        #print('\t Mean     Std') #disp ([' Mean'  '     Std']);
        
        formatSpec = '{:4.1f}' #formatSpec = '%4.1f';
        for i in range(len(sel2)): #disp ([num2str(100*mean(sel2_solutions,2),formatSpec) percent plusMinus6 num2str(100*stDev,formatSpec) percent;]);
            strFormat = SStruct[i]+'\t'+ formatSpec+' +/-' + formatSpec
            writeStr = strFormat.format(100*sel2[i], 100*stDev[i])
            return_List.append(writeStr)
            #print(writeStr)
            #print(strFormat.format(100*sel2[i], 100*stDev[i])) 
        writeStr = 'sum =\t' + formatSpec.format(100*np.sum(sel2)) + '%'
        return_List.append(writeStr)
        #print(writeStr)
        #print('sum =\t', formatSpec.format(100*np.sum(sel2)), '%')   #disp(['sum = ' num2str(sum(100*mean(sel2_solutions,2)),formatSpec) '%']);
        
        SStruct_for_df = SStruct + ['Total']
        sel2_for_df    = np.append(sel2, np.sum(sel2))

        results_df_selcon = pd.DataFrame({'Component_Selcon2':SStruct_for_df,'Percentage':np.round(sel2_for_df*100,1)})

    #end
    
    #%Now refit Selcon solutions, SVH 25/9-06
    
    #print('refit_np_param ', refit_np_param)
    refit_prot = np.zeros((A.shape[0],len(refit_np_param))) #number of datapoints (66) x number of np for refit
    for i in range(len(refit_np_param)): #for(i=1:size(refit_np_param,2))
        #Note the values stores in np_params are the i and j indices
        np_=refit_np_param[i][0] #np_=refit_np_param(1,i);
        basis=refit_np_param[i][1] #basis=refit_np_param(2,i);
        #Note: np is in the range 0 to 71-1, and basis is in the range 0 to 7-1. So use them as indices as they are (not +1)
        u, s, vh = np.linalg.svd(A[:,:np_+1], full_matrices=True) #[U,S,V] = svd(A(:,1:np),0);         #Vt=V';
        S = np.zeros((u.shape[0],vh.shape[0])) 
        S[:s.shape[0],:s.shape[0]] = np.diag(s)
        refit = u[:,:basis+1] @ S[:basis+1,:basis+1] @ vh[:basis+1,:] #        refit = U(:,1:basis)*S(1:basis,1:basis)*Vt(1:basis,:);
        refit_prot[:,i]=refit[:,0] #        refit_prot(:,i)=refit(:,1); %This is the refit prot spec for given #Np and #Basis           

    #end
    #print(refit_prot) 

    Mean_refit_prot = np.mean(refit_prot,axis=1) #mean(refit_prot,2); %Vector of mean refitted prot spectrum

    #pred.MeanRefitProt=Mean_refit_prot;
    #pred.RefitProt=refit_prot; %Array of all refited spectra for all good Sel3 #Np and #Basis
    
    #%Calc. RMSD value
    rmsdRefit = np.sqrt(((Mean_refit_prot - q) ** 2).mean()) #rmsdRefit=rmsd(Mean_refit_prot,q);
    writeStr = 'RMSD between Query and Refitted spectra ' + '{:.4f}'.format(rmsdRefit)
    return_List.append(writeStr)
    ##print(writeStr)
    #print('RMSD between Query and Refitted spectra ', '{:.4f}'.format(rmsdRefit)) #disp(['RMSD between Query and Refitted spectra ' num2str(rmsdRefit)]);
    
    #np.savetxt('Mean_refit_prot.txt', Mean_refit_prot)
    
    NpMax = 0
    
    for i in range(refit_prot.shape[1]):
        if (refit_np_param[i][0]+1 > NpMax):
            NpMax = refit_np_param[i][0]+1
    
    if False:

        #%Calc. selcon2 solutions
        refit_prot_Sel2 = np.zeros((A.shape[0],len(np_param_Sel2))) #number of datapoints (66) x number of np for refit
        for i in range(len(np_param_Sel2)): #for(i=1:size(refit_np_param,2))
            np_   = np_param_Sel2[i][0] #np=np_param_Sel2(1,i);
            basis = np_param_Sel2[i][1] #basis=np_param_Sel2(2,i);
            #Note: np is in the range 0 to 71-1, and basis is in the range 0 to 7-1. So use them as indices as they are (not +1)
            u, s, vh = np.linalg.svd(A[:,:np_+1], full_matrices=True) #[U,S,V] = svd(A(:,1:np),0);         #Vt=V';
            S = np.zeros((u.shape[0],vh.shape[0])) 
            S[:s.shape[0],:s.shape[0]] = np.diag(s)
            refit = u[:,:basis+1] @ S[:basis+1,:basis+1] @ vh[:basis+1,:] #        refit = U(:,1:basis)*S(1:basis,1:basis)*Vt(1:basis,:);
            refit_prot_Sel2[:,i]=refit[:,0] #        refit_prot(:,i)=refit(:,1); %This is the refit prot spec for given #Np and #Basis  
        
    #end %if(ShouldIPlot==1) %So plot solutions
    
    #Which proteins where included in the solutions
    #print('')
    #print('\n****** Which proteins used for the calculation? ******')

    writeStr ='\n'
    writeStr +='****** Which proteins used for the calculation? ******'+'\n'
    writeStr += 'Max Np ' + '{:d}'.format(NpMax)+'\n'
    
    try:
            
        writeStr += 'Proteins used in solution in RMSD order:\n'
        for i in range(NpMax):
            writeStr += '{:d}: '.format(ix[i])
            writeStr += Lbl[ix[i]]+'\n'
        #print(writeStr)
    except:
        writeStr +="Could not find the Label file"
    #END: Try/Except block

    if (includeProteinsUsedForCalc>=1):
        return_List.append(writeStr)
    
    return return_List, results_df_selcon, Mean_refit_prot, query_spectrum
        