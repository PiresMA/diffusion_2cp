/** Authors: M.A.Pires and Queiros
*** Main program to the Monte Carlo Simulation of the model presented in paper
*** Optimal diffusion in ecological dynamics  with Allee effect in a metapopulation
***
*** Specifically this code yields the results shown in Fig.3-a, but with small
*** modifications all the other results presented in our paper can be obtained.
*** Do not hesitate to contact me: piresma@cbpf.br  or pires.ma.fisica@gmail.com
***/

// NICE: http://wakingup.libsyn.com/145-the-information-war

#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ran2_new.h"


int simu=4;

#define nsubp 10
#define sizesubpop 100000
#define N sizesubpop*nsubp    //int N=1000*nsubp;
#define L nsubp               //int nn=2;
#define nn 2

#define tmax      201 //number of time steps to simulate over
#define t_steady  1   // NAO USO AQUI: time to reach the steady state
#define cnt       2   // distancia entre medidas

#define lam 1.0       //
#define r 0.86        //

#define fao 0.48
#define fbo 0.52

int state[N], state2[N];
int local[N], local2[N];





/***Choosing a second agent Inside Subpopulation ***/
int secondAgentInsideSubpopulation(int i1, long *Seed)
{
    int flag=0, i2, possible_i2=i1;

    //if( local_population[i]>2 ) as the S and I do not move all the subpopulation have more than 2 agents
    while(flag==0)
    {
        possible_i2 = ((int)N*ran2(Seed));
        if(possible_i2!=i1)
        {
            if( local[possible_i2]==local[i1] )
            {
                flag=1;
                i2=possible_i2;
                //printf("i:%d site:%d i2:%d site2:%d     ",i1, state[possible_i1],i2,state[possible_i2]);
                return(i2);
            }
        }
    }
}





//int main(int argc, char *argv[])
int main(int argc, char *argv[])
{

    /**********************************************************/
    long SEED=time(NULL);
    if (SEED > 0) SEED = -SEED; /* Seed passed to ran2() must initially be negative. */
    //Seed: a negative number means start a new sequence of
    //pseudo-random numbers; a positive number means continue
    //with the same sequence.  S is turned positive by ran2().
    /**********************************************************/


    int pop,i,j,jj,label,count,soma=0;
    int subpop,index,newsubpop;
    int xx,aux,cont,cont2,tot;
    int site,agent,t;
    int node1,node2,node3;
    int A_cont[L+1],V_cont[L+1],B_cont[L+1];
    int vecn[L],vecAo[L],vecVo[L],vecBo[L];
    int neig_mat[L][nn+1];



    /*****************************************************************************/
    /*****************************************************************************/
    // get int x and convert it into real/float
    // double D = 0.0; // #mobility parameter
    // int arg1=atof(argv[1]);

    //double vecD[9]; // $D=\{0.03,0.05,0.07,0.09,0.14,0.19\}$

    double D   =  0.1*atof(argv[1]);
    double alp =  0.1*atof(argv[2]);
    int nSourc =    1*atof(argv[3]);
    int sample =    1*atof(argv[4]);

    //ktot = 2; // 2*atof(argv[2]);
    /*****************************************************************************/
    /*****************************************************************************/


    /************************************************************************/
    int Aoglobal=0, Voglobal=0, Boglobal=0;
    int auxA2, auxB2;
    auxA2 = (int) ( (fao)*sizesubpop/nSourc );
    auxB2 = (int) ( (fbo)*sizesubpop/nSourc);
    printf("N:%d   auxA2:%d    auxB2:%d \n", N, auxA2, auxB2);
    for(jj=0; jj<nSourc; jj++)
    {
        vecn[jj]  = sizesubpop;
        vecAo[jj] = auxA2;  // maxIo*(L-0)=(1/L)*L=1    //aux2*fatorIo*maxIo*(L-0);
        vecBo[jj] = auxB2;  // maxIo*(L-0)=(1/L)*L=1    //aux2*fatorIo*maxIo*(L-0);
        vecVo[jj] = vecn[jj] - vecAo[jj] - vecBo[jj];
        Aoglobal += 1.0*vecAo[jj];
        Boglobal += 1.0*vecBo[jj];
        Voglobal += 1.0*vecVo[jj];
    }
    for(jj=nSourc; jj<L; jj++)
    {
        vecn[jj]  = sizesubpop;
        vecAo[jj] = 0;
        vecBo[jj] = 0;
        vecVo[jj] = vecn[jj] - vecAo[jj] - vecBo[jj];
        Aoglobal += 1.0*vecAo[jj];
        Boglobal += 1.0*vecBo[jj];
        Voglobal += 1.0*vecVo[jj];
    }
    printf("Aoglobal:%d   Voglobal:%d  Boglobal:%d \n", Aoglobal, Voglobal, Boglobal);
    printf("Aoglobal:%1.3lf  Voglobal:%1.3lf  Boglobal:%1.3lf\n", Aoglobal/(1.0*N), Voglobal/(1.0*N), Boglobal/(1.0*N) );

    /************************************************************************/


    printf("neig_mat:\n");
    for(i=0; i<L; i++)// condicoes de contorno periodicas em x:
    {
        neig_mat[i][0] = i;
        neig_mat[i][1] = (i+L-1)%L; //im=(i+L-1)%L;
        neig_mat[i][2] = (i+1)%L; //ip=(i+1)%L;
        printf("%d  %d  %d\n",neig_mat[i][0],neig_mat[i][1],neig_mat[i][2]);
    }



    /**********************************************/
    /**********************************************/
    // make data directory
    char dir[500];
    sprintf(dir,"2cp-simu%d",simu);
    mkdir(dir, 0755);// 0755 is related to Folder/Directory Permissions.
    //askubuntu.com/questions/638796/what-is-meaning-of-755-permissions-in-samba-share

    FILE *fM;
    char outNameM[500];
    /**********************************************/
    /**********************************************/

    //  int nSourc =    1*atof(argv[3]);
    sprintf(outNameM,"./%s/2cp_simu%d_nS%d_r%1.3lf_alp%1.3lf_N%d_L%d_D%1.2lf_tmax%d_samp%d.dat",dir,simu,nSourc,r,alp,N,L,D,tmax,sample);
    fM=fopen(outNameM,"w");

    /************ Initialization   ***********/

    label=0;
    for(jj=0; jj<L; jj++)
    {
        if(vecAo[jj]!=0)
        {
            for(i=0; i<vecAo[jj]; i++)
            {
                state[label] = -1;                
                local[label] = jj;

                state2[label]= state[label];
                local2[label]= local[label];

                label++;
            }
        }
        
        if(vecVo[jj]!=0)
        {
            for(i=0; i<vecVo[jj]; i++)
            {
                state[label] = 0;                
                local[label] = jj;

                state2[label]= state[label];
                local2[label]= local[label];

                label++;
            }
        }

        if(vecBo[jj]!=0)
        {
            for(i=0; i<vecBo[jj]; i++)
            {
                state[label] = 1;                
                local[label] = jj;

                state2[label]= state[label];
                local2[label]= local[label];

                label++;
            }
        }
    }

        for(pop=0; pop<L; pop++)
        {
            A_cont[pop]=0;
            V_cont[pop]=0;
            B_cont[pop]=0;
        }

        for(i=0; i<N; i++)
        {
                if( state[i] == -1 ) A_cont[local[i]]++;
                else if( state[i] == 0 ) V_cont[local[i]]++;
                else B_cont[local[i]]++;
        }



/****
    printf("state:\n");
    for(j=0; j<N; j++) printf("%d  ",state[j]);
    printf(" \n");

    printf("local:\n");
    for(j=0; j<N; j++) printf("%d  ",local[j]);
    printf(" \n");


    for(pop=0; pop<L; pop++) printf("t:%d A_cont[%d]:%d V_cont[%d]:%d B_cont[%d]:%d\n",t,pop,A_cont[pop],pop,V_cont[pop],pop,B_cont[pop]);
    printf("\n");
***/


    fprintf(fM,"%d  ",0);
    for(pop=0; pop<L; pop++) fprintf(fM,"%d %d %d   ",A_cont[pop],V_cont[pop],B_cont[pop]);
    fprintf(fM,"\n");

    /************ PRINT AND CREATE NAMEFILE  ***********/



    printf("****SEED:%ld  fao:%1.2lf  fbo:%1.2lf  r:%1.3lf  lam:%1.1lf alp:%1.2lf  D:%1.2lf  L:%d  ttot:%d\n",SEED,fao,fbo,r,lam,alp,D,L,tmax);
    //printf("*** SEED:%ld  simu:%d  L:%d  ttot:%d\n",SEED,simu,L,tmax);

    /************************   Dynamics (1): BEGIN   ***********************/
    count  = 0;
    for(t=1; t<tmax; t++)
    {
        for(node1=0; node1<N; node1++)
        {
            if(  ran2(&SEED) < D ) //ran2(Seed)<D
            {
                if(  state[node1]!=0 ) // if i is A or B
                {
                    subpop        = local[node1];
                    index         =  1 + (int) ( ran2(&SEED)*nn );  //A+(int)(ran2(Seed)*(B-A))-->[A,B-1]//range [1,nn]
                    newsubpop     = neig_mat[subpop][index];
                    local2[node1] = newsubpop;
                    //printf("t:%d %d  from:%d  index:%d  to:%d  %d\n",t,node1,subpop,index,newsubpop,state[node1]);
                }
            }
            else
            {
                if(  state[node1] == 0 ) // if i is V
                {
                   soma=A_cont[local[node1]]+B_cont[local[node1]];
                    if( soma>0 ) node2=secondAgentInsideSubpopulation(node1,&SEED);
                    else node2=node1;
                                        
                    
                    if( state[node2] == -1 ) // if j is A
                    {
                        if( ran2(&SEED) < lam ) state2[node1] = -1; // V+A-->2A                   
                    }
                    else if( state[node2] == 1 ) // if j is B
                    {
                        if( ran2(&SEED) < r*lam ) state2[node1] = 1; // V+B-->2B
                    }
                   //printf("nodes: %d  %d  states:%d  %d  soma:%d \n",node1,node2,state[node1],state[node2],soma);
                }   
                else
                {
                    if( state[node1] == -1 ) // if i is A
                    {
                        if( ran2(&SEED) < alp ) state2[node1]=0; // A-->V
                    }
                    else if( state[node1] == 1 ) // if i is B
                    {
                        if( ran2(&SEED) < r*alp ) state2[node1]=0; // B-->V
                    }
                }
            }

        }//End of MC step over all the N states

        for(pop=0; pop<L; pop++)
        {
            A_cont[pop]=0;
            V_cont[pop]=0;
            B_cont[pop]=0;
        }

        for(i=0; i<N; i++)
        {
            state[i]=state2[i]; // Parallel/Syncronous Updating
            local[i]=local2[i]; // Parallel/Syncronous Updating

                if( state[i] == -1 ) A_cont[local[i]]++;
                else if( state[i] == 0 ) V_cont[local[i]]++;
                else B_cont[local[i]]++;
        }

 

    //for(pop=0; pop<L; pop++) printf("t:%d A_cont[%d]:%d V_cont[%d]:%d B_cont[%d]:%d\n",t,pop,A_cont[pop],pop,V_cont[pop],pop,B_cont[pop]);
    //printf("\n");


        count++;
        if(count==cnt)
        {
            count=0;
            fprintf(fM,"%d  ",t);
            for(pop=0; pop<L; pop++) fprintf(fM,"%d %d %d   ",A_cont[pop],V_cont[pop],B_cont[pop]);
            fprintf(fM,"\n");
        }

       
    }// END time evolution
    /************************   Dynamics: END   ***********************/


    printf("END!!!!!!!!!!!!\n");




    return(0);
}


