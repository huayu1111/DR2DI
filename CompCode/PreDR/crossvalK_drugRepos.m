function [prdY,Yr]  = crossvalK_drugRepos( Adj, sim_drug, sim_disease, fold )

% the row of Adj represents drug, the column of Adj represents disease;
% sim_drug is drug similarity matrix; sim_disease is disease similarity
% prdY is the probility of label which is used to plot ROC curve;
% Yr is the primal output with the value of {1,-1};

Yr=[]; prdY = [];

[Ip,Jp]=find(Adj==1);
[Iu,Ju]=find(Adj==0);
p=randperm(length(Iu));
In=Iu(p(1:length(Ip)));Jn=Ju(p(1:length(Ip)));

mp=length(Ip);
p1=randperm(mp);
Ipr=Ip(p1);Jpr=Jp(p1);

mn=length(In);
p2=randperm(mn);
Inr=In(p2);Jnr=Jn(p2);

mpr=length(Ipr);
mnr=length(Inr);

t = fix( mpr / fold);

for k = 1 : fold
    Iprn=[];Jprn=[];
    Iurn=[];Jurn=[];
    Ipst=[];Jpst=[];
    Inst=[];Jnst=[];
    if k==fold
        Iprn = Ipr(1 : t * ( k - 1 ));
        Jprn = Jpr(1 : t * ( k - 1 ));
        Inrn = Inr(1 : t * ( k - 1 ));
        Jnrn = Jnr(1 : t * ( k - 1 ));
        Ipst = Ipr( t * ( k - 1 ) + 1 : mpr);
        Jpst = Jpr( t * ( k - 1 ) + 1 : mpr);
        Inst = Inr( t * ( k - 1 ) + 1 : mnr); 
        Jnst = Jnr( t * ( k - 1 ) + 1 : mnr);
   else
       Iprn = [ Ipr( 1 : t* ( k - 1 ) ); Ipr( t * k + 1 : mpr ) ];
       Jprn = [ Jpr( 1 : t * ( k - 1 ) ); Jpr( t * k + 1 : mpr ) ];
       Inrn = [ Inr( 1 : t * ( k - 1 ) ); Inr( t * k + 1 : mnr ) ];
       Jnrn = [ Jnr( 1 : t * ( k - 1 ) ); Jnr( t * k + 1 : mnr ) ];
       Ipst = Ipr( t * ( k - 1 ) + 1 : t * k );
       Jpst = Jpr( t * ( k - 1 ) + 1 : t * k );
       Inst = Inr( t * ( k - 1 ) + 1 : t * k );
       Jnst = Jnr( t * ( k - 1 ) + 1 : t * k );
    end   
    LabGrn=[Iprn;Inrn];LabDrn=[Jprn;Jnrn];
    LabGst=[Ipst;Inst];LabDst=[Jpst;Jnst];
    for i=1:length(LabGrn)
        Yrn(i)=Adj(LabGrn(i),LabDrn(i));
    end
    Y2=Yrn';Y2(Y2==0)=-1;
    for i=1:length(LabGrn)
        for j=1:length(LabDrn)
            Krn(i,j)=sim_drug(LabGrn(i),LabGrn(j))*sim_disease(LabDrn(i),LabDrn(j));
        end
    end
    Drn=[Y2 Krn];   formatKmatrix(Drn,'Dtrain.txt');
    
    for i=1:length(LabGst)
        Yst(i)=Adj(LabGst(i),LabDst(i));
    end
    Y3=Yst';Y3(Y3==0)=-1;Yr=[Yr;Y3];
    for i=1:length(LabGst)
        for j=1:length(LabDrn)
            Kst(i,j)=sim_drug(LabGst(i),LabGrn(j))*sim_disease(LabDst(i),LabDrn(j));
        end
    end
%     [Hst] = Normal_M(Kst);
    Dst=[Y3 Kst];   formatKmatrix(Dst,'Dtest.txt');
      
  !E:\drug2disease\DR2DI\Data\建模数据\CompCurve\PreDR\svm-train.exe -s 0 -b 1 -t 4 -c 1.4 Dtrain.txt 
  % make sure software svm-train is under this directory
      
  !E:\drug2disease\DR2DI\Data\建模数据\CompCurve\PreDR\svm-predict.exe -b 1 Dtest.txt Dtrain.txt.model result.out
  % make sure software svm-predict is under this directory
  
   [names,types,x] = textread('result.out','%s%s%s');
   if strcmp( types{1} , '-1' )
       types = x;
   end
   
   predictedY = zeros( length(types) - 1 , 1 );
   for i=2:length(types)
       predictedY(i-1)=str2double(types{i});
   end
   prdY=[prdY;predictedY];    
end
