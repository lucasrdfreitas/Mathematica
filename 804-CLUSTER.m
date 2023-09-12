(* ::Package:: *)

t0=AbsoluteTime[];
Print["Before Starting Kernels"];


Needs["ClusterIntegration`"];
(*kernels = aunchKernels[SGE["micro4", 10]];*)
Quiet[kernels = LaunchKernels[]];
Print["Starting Kernels"];


(* ::Section::Bold:: *)
(*Definitions *)


NbName="804"; \[Lambda]0=0.5; 

Ls =Range[40,40,2]; 	    
tV={0};	  
hV={{.2,0,0}};
steps=150;
acuracy=7;    
\[CapitalDelta]eV=0.1; eVs=Table[ 1700 \[Xi], {\[Xi],0,0,\[CapitalDelta]eV} ] ;
\[Eta]s=Join[ Table[3 \[Eta]-.113,{\[Eta],-1,1,.1}] (*, Table[\[Eta],{\[Eta],-.5+.02,.5,0.05}], Table[\[Eta],{\[Eta],-.5+.04,.5,0.05}] *) ];
		


(* ::Subsubsection::Bold::Closed:: *)
(*Constants*)


ReverseSort[list_]:=Reverse@Sort@list;



one[i_,j_,N_]:= SparseArray[ {i,j} -> 1, {N,N}];  block[i_,j_] := one[i,j,8];
hc=1/Sqrt[3] {1,1,1};hb=1/Sqrt[2] {1,-1,0};ha=1/Sqrt[6] {1,1,-2};hx={1,0,0};
 hAngle[\[Theta]_,\[Phi]_]:= ha Cos[\[Phi] \[Pi]/180] Sin[\[Theta] \[Pi]/180]+hb Sin[\[Phi] \[Pi]/180] Sin[\[Theta] \[Pi]/180]+ hc Cos[\[Theta] \[Pi]/180] ;  (* Magnetic field directions a, b, c*)
nx={1/2,Sqrt[3]/2};  ny={-(1/2),Sqrt[3]/2};                                                                                           (* Honeycomb Bravais basis  *)

\[Delta]x={-Sqrt[3]/2,-1/2}/Sqrt[3];      \[Delta]y={Sqrt[3]/2,-1/2}/Sqrt[3];        \[Delta]z={0,1}/Sqrt[3];   

asites[m_,n_]:=m nx+n ny;
bsites[m_,n_]:=m nx+n ny-\[Delta]z;
\[Sigma]sites[m_,n_,\[Sigma]_]:=m nx+n ny-\[Sigma] \[Delta]z;
Gmat=  {KroneckerProduct[  PauliMatrix[0],-I PauliMatrix[2] ],KroneckerProduct[- I PauliMatrix[2],PauliMatrix[3] ],KroneckerProduct[-I PauliMatrix[2],  PauliMatrix[1] ]};
Mmat=  {KroneckerProduct[  PauliMatrix[3],I PauliMatrix[2] ],KroneckerProduct[ I PauliMatrix[2],PauliMatrix[0] ],KroneckerProduct[  PauliMatrix[1],I PauliMatrix[2] ]};
Nmat=(Mmat-Gmat)/2;
to\[Lambda][h_,\[CapitalDelta]v_:0.262]:= {h[[1]]^2,h[[2]]^2,h[[3]]^2}/\[CapitalDelta]v;
toKappa[h_,\[CapitalDelta]v_:0.262]:=8h[[1]]h[[2]]h[[3]]/(  3 \[CapitalDelta]v^2 );
KappaToH[\[Kappa]_,d_,\[CapitalDelta]v_:0.262]:=Module[{C=d[[1]]d[[2]]d[[3]]},If[C==0,{0,0,0},  d CubeRoot[3  \[CapitalDelta]v^2 \[Kappa]/(8C)]   ]] ;
\[Chi]Gx = {{0.5249,0,0,0},{0,-1,0,0},{0,0,0,0},{0,0,0,0}};\[Chi]Gy = {{0.5249,0,0,0},{0,0,0,0},{0,0,-1,0},{0,0,0,0}};\[Chi]Gz ={{0.5249,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,-1}};
\[Omega]GA = I IdentityMatrix[4]+10^-6 Sum[Mmat[[\[Gamma]]],{\[Gamma],1,3}]; 
\[Omega]GB = \[Omega]GA;
Uguess=-{\[Chi]Gx,\[Chi]Gy,\[Chi]Gz};
Vguess={\[Omega]GA,\[Omega]GB};
to3digits[list0_]:=Module[{list=list0},If[ Length@list>=3, list, list=Flatten[{{0},list},1]; If[ Length@list>=3, list,  list=Flatten[{{0},list},1] ] ]; StringJoin[ToString/@list]     ];
round[x_]:=N[Round[1000000(x)]/1000000];
(* Module[{V=Table[v[[i,j]],{i,1,4},{j,1,4}]  },Table[Tr[V^\[Transpose].Nmat[[\[Gamma]]]],{\[Gamma],1,3}]] *)
sumtraceG[v_]:=-v[[1,2]]-v[[1,3]]-v[[1,4]]+v[[2,1]]-v[[2,3]]+v[[2,4]]+v[[3,1]]+v[[3,2]]-v[[3,4]]+v[[4,1]]-v[[4,2]]+v[[4,3]];
traceG[v_]:={-v[[1,2]]+v[[2,1]]-v[[3,4]]+v[[4,3]],-v[[1,3]]+v[[2,4]]+v[[3,1]]-v[[4,2]],-v[[1,4]]-v[[2,3]]+v[[3,2]]+v[[4,1]]};
sumtraceM[v_]:=v[[1,2]]-v[[2,1]]-v[[3,4]]+v[[4,3]]+v[[1,3]]+v[[2,4]]-v[[3,1]]-v[[4,2]]+v[[1,4]]-v[[2,3]]+v[[3,2]]-v[[4,1]];
traceM[v_]:={v[[1,2]]-v[[2,1]]-v[[3,4]]+v[[4,3]],v[[1,3]]+v[[2,4]]-v[[3,1]]-v[[4,2]],v[[1,4]]-v[[2,3]]+v[[3,2]]-v[[4,1]]};
sumtraceN[v_]:=v[[1,2]]+v[[1,3]]+v[[1,4]]-v[[2,1]]-v[[3,1]]-v[[4,1]];
traceN[v_]:={v[[1,2]]-v[[2,1]],v[[1,3]]-v[[3,1]],v[[1,4]]-v[[4,1]]};
traceUNUN[u_]:=2{{u[[1,2]] u[[2,1]]-u[[1,1]] u[[2,2]],u[[1,3]] u[[2,1]]-u[[1,1]] u[[2,3]],u[[1,4]] u[[2,1]]-u[[1,1]] u[[2,4]]},
{u[[1,2]] u[[3,1]]-u[[1,1]] u[[3,2]],u[[1,3]] u[[3,1]]-u[[1,1]] u[[3,3]],u[[1,4]] u[[3,1]]-u[[1,1]] u[[3,4]]},{u[[1,2]] u[[4,1]]-u[[1,1]] u[[4,2]],u[[1,3]] u[[4,1]]-u[[1,1]] u[[4,3]],u[[1,4]] u[[4,1]]-u[[1,1]] u[[4,4]]}};



(* ::Subsection::Bold::Closed:: *)
(*for pure*)


to\[Lambda][h_,\[CapitalDelta]v_:0.262]:= {h[[1]]^2,h[[2]]^2,h[[3]]^2}/(2\[CapitalDelta]v);
toKappa[h_,\[CapitalDelta]v_:0.262]:=8h[[1]]h[[2]]h[[3]]/(  3 \[CapitalDelta]v^2 );
KappaToH[\[Kappa]_,d_,\[CapitalDelta]v_:0.262]:=Module[{C=d[[1]]d[[2]]d[[3]]},If[C==0,{0,0,0},  d CubeRoot[3  \[CapitalDelta]v^2 \[Kappa]/(8C)]   ]] ;


(* ::Subsubsection::Bold::Closed:: *)
(*file*)


dataToFilePure[ parameters0_,L_,acuracy_,gauge_,data_] :=
Module[ {path,f,parameters=parameters0},
parameters[[1]]=0;parameters[[3]]=0;parameters[[5]]=0;parameters[[7]]=0;		
(*createDir@FileNameJoin[{Directory[],"Files","pure", gauge}] ;*)
createDir@FileNameJoin[{Directory[],"Files",NbName,"pure",gauge, StringReplace["t=X1_eV=X2_JKG=X3_JKGmod=X4",
{"X1"->ToString[parameters[[9]]],"X2"->ToString[parameters[[10]]],"X3"->ToString[parameters[[1;;3]]],"X4"->ToString[parameters[[5;;7]]] }] }];
		path = toPathPure[parameters,L,acuracy,gauge,NbName];		
		(*Print["Pure path=",path];Print[];*)
		f = OpenAppend[path];
		 Write[ f, data];
		 Close[f];                ];	
		 		 
toPathPure[parameters0_,L_,acuracy_,gauge_]:= 
 Module[{h ,hS,parameters=parameters0, r,\[Phi],\[Theta]},h=parameters[[4]] ;hS=parameters[[8]]; {r,\[Theta],\[Phi]}=hS;   
parameters[[1]]=0;parameters[[3]]=0;parameters[[5]]=0;parameters[[7]]=0;
FileNameJoin[{Directory[], "Files" ,NbName,"pure",gauge,StringReplace["t=X1_eV=X2_JKG=X3_JKGmod=X4",
{"X1"->  ToString[parameters[[9]]],"X2"->  ToString[parameters[[10]]],"X3"->  ToString[parameters[[1;;3]]  ],
"X4"->  ToString[parameters[[5;;7]]]   }] , "data"  , 
StringReplace["h=(M,N,T)_L=Y_A=Z.txt",{"Y"-> ToString[L], "Z"-> ToString[acuracy],"M"->  ToString[r,InputForm] ,"N"->ToString@\[Phi],"T"->ToString@\[Theta]}   ] 
   }]];
 


createDir[path_] :=
Module[ {l=Length@FileNames[path]},		
	If[ l==0,  (*CreateDirectory@File@FileNameJoin[{path}];*) 
	CreateDirectory@File@FileNameJoin[{path,"data" }]; CreateDirectory@File@FileNameJoin[{path,"graph" }]; , Null]  ];


dataToFilePure[ parameters0_,L_,acuracy_,data_,gauge_,NbName_]:=
Module[ {path,f,parameters=parameters0},
parameters[[1]]=0;parameters[[3]]=0;parameters[[5]]=0;parameters[[7]]=0;		
(*createDir@FileNameJoin[{Directory[],"Files","pure", gauge}] ;*)
createDir@FileNameJoin[{Directory[],"Files",NbName,"pure",gauge, StringReplace["t=X1_eV=X2_JKG=X3_JKGmod=X4",
{"X1"->ToString[parameters[[9]]],"X2"->ToString[parameters[[10]]],"X3"->ToString[parameters[[1;;3]]],"X4"->ToString[parameters[[5;;7]]] }] }];
		path = toPathPure[parameters,L,acuracy,gauge,NbName];		
		(*Print["Pure path=",path];Print[];*)
		f = OpenAppend[path];
		 Write[ f, data];
		 Close[f];                ];		 
		 
toPathPure[parameters0_,L_,acuracy_,gauge_,NbName_]:= Module[{h ,hS,parameters=parameters0, r,\[Phi],\[Theta]},h=parameters[[4]] ;hS=parameters[[8]]; {r,\[Theta],\[Phi]}=hS;   
parameters[[1]]=0;parameters[[3]]=0;parameters[[5]]=0;parameters[[7]]=0;
FileNameJoin[{Directory[], "Files" ,NbName,"pure",gauge,StringReplace["t=X1_eV=X2_JKG=X3_JKGmod=X4",
{"X1"->  ToString[parameters[[9]]],"X2"->  ToString[parameters[[10]]],"X3"->  ToString[parameters[[1;;3]]  ],
"X4"->  ToString[parameters[[5;;7]]]   }] , "data"  , 
StringReplace["h=(M,N,T)_L=Y_A=Z.txt",{"Y"-> ToString[L], "Z"-> ToString[acuracy],"M"->  ToString[r,InputForm] ,"N"->ToString@\[Phi],"T"->ToString@\[Theta]}   ] 
   }]];

loadDataPure[path_]:=Module[ {f,data},
f = OpenRead[path];
If[f==$Failed, Print["Failed to OpenRead file at: "]; Print[ path ]; Abort[] ];
data=ReadList[f];
Close[f];			
data[[-1]]
];



(* ::Subsubsection::Bold::Closed:: *)
(*Auxiliary matrices for the Hamiltonian  [2x2 matrices]*)


(* ::Text::Bold:: *)
(*we define Z2 gauge field u[r] for each unit cell r as a three component vector which each component correspond \[Alpha] to the gauge field in the \[Alpha] bond for the A site in the unit cell r.*)


(* ::Code::Bold:: *)
(**)


toR[m0_,n0_,L_,M_]:= Module[{\[CapitalDelta]n=\[LeftFloor]n0/L\[RightFloor],n,m},n=Mod[n0 ,L]; m=Mod[m0+M \[CapitalDelta]n,L]; m + n L+1];
HNN[K_,\[Kappa]_,\[Lambda]_,u_,m_,n_,L_,M_]:=Module[{r=toR[m,n ,L,M]},Table[{{0,(-K[[r,\[Alpha]]]+\[Lambda][[\[Alpha]]]) u[[r,\[Alpha]]] },{0,0}},{\[Alpha],1,3}]]; 
HNNNA[K_,\[Kappa]_,\[Lambda]_,u_,m_,n_,L_,M_]:=Module[{
r1={toR[m,n ,L,M],toR[m,n ,L,M],toR[m,n,L,M]},
r2={toR[m+1,n-1 ,L,M],toR[m,n+1 ,L,M],toR[m-1,n,L,M]}        }      ,
Table[{{\[Kappa] u[[r1[[\[Alpha]]],\[Alpha]]] u[[r2[[\[Alpha]]],Mod[\[Alpha]+1,3,1] ]],0},{0,0 }},{\[Alpha],1,3}]
];
HNNNB[K_,\[Kappa]_,\[Lambda]_,u_,m_,n_,L_,M_]:=Module[{
r1={toR[m-1,n ,L,M],toR[m,n -1,L,M],toR[m,n,L,M]},
r2={toR[m-1,n ,L,M],toR[m,n-1 ,L,M],toR[m,n,L,M]}        }      ,
Table[{{0,0},{0,\[Kappa] u[[r1[[\[Alpha]]],\[Alpha]]] u[[r2[[\[Alpha]]],Mod[\[Alpha]+1,3,1] ]]}},{\[Alpha],1,3}]
];
					
Hreal[K_,\[Kappa]_,\[Lambda]_,u_,L_,M_,\[Delta]k_] := Module[{Nc=L^2,H },

 H= I Sum[Module[{Hnn,HnnnA,HnnnB, \[CurlyPhi]=\[Delta]k . (m nx+n ny),
r={toR[m+1,n ,L,M],toR[m,n+1 ,L,M],toR[m,n ,L,M]},
RA={toR[m+1,n-1 ,L,M],toR[m,n+1 ,L,M],toR[m-1,n,L,M]},
RB={toR[m-1,n+1 ,L,M],toR[m,n-1 ,L,M],toR[m+1,n,L,M]}    },

Hnn=HNN[K,\[Kappa],\[Lambda],u,m,n,L,M] ;HnnnA=HNNNA[K,\[Kappa],\[Lambda],u,m,n,L,M] ;HnnnB=HNNNB[K,\[Kappa],\[Lambda],u,m,n,L,M] ;

Sum[       KroneckerProduct[ one[r[[3]],r[[\[Alpha]]],Nc],Hnn[[\[Alpha]]]]   +
KroneckerProduct[ one[r[[3]],RA[[\[Alpha]]],Nc],HnnnA[[\[Alpha]]]]   +
KroneckerProduct[ one[r[[3]],RB[[\[Alpha]]],Nc],HnnnB[[\[Alpha]]]]   ,{\[Alpha],1,3}] Exp[I \[CurlyPhi]]      ] 
,{m,0,L-1},{n,0,L-1}];

 H+H\[ConjugateTranspose]];

TmatPure[L_] :=KroneckerProduct[   IdentityMatrix[L^2],  {{1,1},{I,-I}} ];
UmatPure[H_]:= Module[ {R=Quiet@Eigensystem@N[H]},ReverseSort[R\[Transpose]]\[Transpose][[2]]\[Transpose] ];
EandUPure[H_]:= Module[ {R=Transpose@ReverseSort@Transpose@Quiet@Eigensystem@N[H]},{R[[1]],R[[2]]\[Transpose] }];
\[DoubleStruckCapitalU]occupiedPure[TU_,Nc_]:= Drop[Take[TU\[Transpose],-Nc-1],{2}]\[Transpose];


(* ::Code::Bold:: *)
(**)


correlations[U_,u_,L_,M_] :=Module[  { Nc=L^2,\[DoubleStruckCapitalU],\[DoubleStruckCapitalU]h,icc,SS },
\[DoubleStruckCapitalU]=TmatPure[L] . U;\[DoubleStruckCapitalU]h=\[DoubleStruckCapitalU][[;;,-Nc;;-1]];icc=Chop[I  \[DoubleStruckCapitalU]h . \[DoubleStruckCapitalU]h\[ConjugateTranspose] ];

SS[[3]]=ArrayFlatten[ # ,1]& @ Table[ 
Module[{ rz,Io,J},rz=m+(n-1)L-1; Io=Mod[\[Alpha]+2rz,2L^2,1];J=Mod[\[Beta] + 4 +2rz,2L^2,1];
-u[[m,n]] icc[[J,Io]]    ], {m,1,L},{n,1,L},  {\[Alpha],1,4},{\[Beta],1,4} ];

SS[[1]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,rx,mx,Io,J},mx=Mod[m+1,L,1]; rz=m+(n-1)L-1;rx=mx+(n-1)L-1;Io=Mod[\[Alpha]+2rz,2L^2,1];J=Mod[\[Beta] + 4+2rx,2L^2,1];
icc[[J,Io]]      ]   , {m,1,L},{n,1,L},  {\[Alpha],1,4},{\[Beta],1,4} ];
SS[[2]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,ry,my,ny,Io,J}, 
ny=Mod[n+1,L,1];my=If[ny==1,Mod[m-M,L,1] , m];  ry=my+(ny-1)L-1;rz=m+(n-1)L-1;Io=Mod[\[Alpha]+2rz,2L^2,1];J=Mod[\[Beta] + 4+2ry,2L^2,1];
icc[[J,Io]]    ]    , {m,1,L},{n,1,L},  {\[Alpha],1,4},{\[Beta],1,4} ];
SS
];

(* writting the exp values in the form of the MF parameters:   *)
toMFparametersPure[U_,u_,L_,M_] :=Module[  { Nc=L^2,\[DoubleStruckCapitalU],\[DoubleStruckCapitalU]h,icc,SS=Array[Null,3]  },
\[DoubleStruckCapitalU]=TmatPure[L] . U;\[DoubleStruckCapitalU]h=\[DoubleStruckCapitalU][[;;,-Nc;;-1]];icc=Chop[I  \[DoubleStruckCapitalU]h . \[DoubleStruckCapitalU]h\[ConjugateTranspose] ];

SS[[3]]=ArrayFlatten[ # ,1]& @ Table[ 
Module[{ rz,Io,J},rz=m+n L; Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1 +2rz,2L^2,1];
 {{    1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0.01,0,0},{0,0,0.01,0},{0,0,0,u[[ Mod[rz+1, L^2,1] ,3]] }}    ],{n,0,L-1} , {m,0,L-1}];

SS[[1]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,rx,mx,Io,J},mx=Mod[m+1,L]; rz=m+n L;rx= mx+n L ;Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1+2rx,2L^2,1];
 {{   1/2 (icc[[J,Io]]-icc[[Io,J]] ) ,0,0,0},{0,u[[ Mod[rz+1,L^2,1] ,1]] ,0,0},{0,0,0.01,0},{0,0,0,0.01}}     ] ,{n,0,L-1} , {m,0,L-1}];
SS[[2]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,ry,my,ny,Io,J}, ny=Mod[n+1,L];my=If[ny==0,Mod[m-M,L] , m];  ry= my+ny L;
rz=m+n L;Io=Mod[1+2rz,2L^2,1];J=Mod[1+ 1+2ry,2L^2,1];
 {{     1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0.01,0,0},{0,0,u[[ Mod[rz+1,L^2,1],2]] ,0},{0,0,0,0.01}}        ],{n,0,L-1} , {m,0,L-1}];
SS        ];

toMFparametersOccupiedPure[U_,u_,L_,M_] :=Module[  { Nc=L^2,\[DoubleStruckCapitalU],\[DoubleStruckCapitalU]h,icc,SS=Array[Null,3]  },
\[DoubleStruckCapitalU]=TmatPure[L] . U;(*\[DoubleStruckCapitalU]h=\[DoubleStruckCapitalU][[;;,-Nc;;-1]];*)\[DoubleStruckCapitalU]h=\[DoubleStruckCapitalU]occupiedPure[\[DoubleStruckCapitalU],Nc];icc=Chop[I  \[DoubleStruckCapitalU]h . \[DoubleStruckCapitalU]h\[ConjugateTranspose] ];

SS[[3]]=ArrayFlatten[ # ,1]& @ Table[ 
Module[{ rz,Io,J},rz=m+n L; Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1 +2rz,2L^2,1];
 {{    1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,u[[ Mod[rz+1, L^2,1] ,3]] }}    ],{n,0,L-1} , {m,0,L-1}];

SS[[1]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,rx,mx,Io,J},mx=Mod[m+1,L]; rz=m+n L;rx= mx+n L ;Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1+2rx,2L^2,1];
 {{   1/2 (icc[[J,Io]]-icc[[Io,J]] ) ,0,0,0},{0,u[[ Mod[rz+1,L^2,1] ,1]] ,0,0},{0,0,0,0},{0,0,0,0}}     ] ,{n,0,L-1} , {m,0,L-1}];
SS[[2]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,ry,my,ny,Io,J}, ny=Mod[n+1,L];my=If[ny==0,Mod[m-M,L] , m];  ry= my+ny L;
rz=m+n L;Io=Mod[1+2rz,2L^2,1];J=Mod[1+ 1+2ry,2L^2,1];
 {{     1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0,0,0},{0,0,u[[ Mod[rz+1,L^2,1],2]] ,0},{0,0,0,0}}        ],{n,0,L-1} , {m,0,L-1}];
SS        ];




(* ::Subsubsection::Bold::Closed:: *)
(*to MF pure*)


(* ::Input:: *)
(*(* writting the exp values in the form of the MF parameters:   *)*)
(*toMFparametersPure[U_,u_,L_,M_] :=Module[  { Nc=L^2,\[DoubleStruckCapitalU],\[DoubleStruckCapitalU]h,icc,SS=Array[Null,3]  },*)
(*\[DoubleStruckCapitalU]=TmatPure[L] . U;\[DoubleStruckCapitalU]h=\[DoubleStruckCapitalU][[;;,-Nc;;-1]];icc=Chop[I  \[DoubleStruckCapitalU]h . \[DoubleStruckCapitalU]h\[ConjugateTranspose] ];*)
(**)
(*SS[[3]]=ArrayFlatten[ # ,1]& @ Table[ *)
(*Module[{ rz,Io,J},rz=m+n L; Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1 +2rz,2L^2,1];*)
(* {{    1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,u[[ Mod[rz+1, L^2,1] ,3]] }}    ],{n,0,L-1} , {m,0,L-1}];*)
(**)
(*SS[[1]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,rx,mx,Io,J},mx=Mod[m+1,L]; rz=m+n L;rx= mx+n L ;Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1+2rx,2L^2,1];*)
(* {{   1/2 (icc[[J,Io]]-icc[[Io,J]] ) ,0,0,0},{0,u[[ Mod[rz+1,L^2,1] ,1]] ,0,0},{0,0,0,0},{0,0,0,0}}     ] ,{n,0,L-1} , {m,0,L-1}];*)
(*SS[[2]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,ry,my,ny,Io,J}, ny=Mod[n+1,L];my=If[ny==0,Mod[m-M,L] , m];  ry= my+ny L;*)
(*rz=m+n L;Io=Mod[1+2rz,2L^2,1];J=Mod[1+ 1+2ry,2L^2,1];*)
(* {{     1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0,0,0},{0,0,u[[ Mod[rz+1,L^2,1],2]] ,0},{0,0,0,0}}        ],{n,0,L-1} , {m,0,L-1}];*)
(*SS        ];*)
(**)
(*toMFparametersOccupiedPure[U_,u_,L_,M_] :=Module[  { Nc=L^2,\[DoubleStruckCapitalU],\[DoubleStruckCapitalU]h,icc,SS=Array[Null,3]  },*)
(*\[DoubleStruckCapitalU]=TmatPure[L] . U;(*\[DoubleStruckCapitalU]h=\[DoubleStruckCapitalU][[;;,-Nc;;-1]];*)\[DoubleStruckCapitalU]h=\[DoubleStruckCapitalU]occupiedPure[\[DoubleStruckCapitalU],Nc];icc=Chop[I  \[DoubleStruckCapitalU]h . \[DoubleStruckCapitalU]h\[ConjugateTranspose] ];*)
(**)
(*SS[[3]]=ArrayFlatten[ # ,1]& @ Table[ *)
(*Module[{ rz,Io,J},rz=m+n L; Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1 +2rz,2L^2,1];*)
(* {{    1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,u[[ Mod[rz+1, L^2,1] ,3]] }}    ],{n,0,L-1} , {m,0,L-1}];*)
(**)
(*SS[[1]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,rx,mx,Io,J},mx=Mod[m+1,L]; rz=m+n L;rx= mx+n L ;Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1+2rx,2L^2,1];*)
(* {{   1/2 (icc[[J,Io]]-icc[[Io,J]] ) ,0,0,0},{0,u[[ Mod[rz+1,L^2,1] ,1]] ,0,0},{0,0,0,0},{0,0,0,0}}     ] ,{n,0,L-1} , {m,0,L-1}];*)
(*SS[[2]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,ry,my,ny,Io,J}, ny=Mod[n+1,L];my=If[ny==0,Mod[m-M,L] , m];  ry= my+ny L;*)
(*rz=m+n L;Io=Mod[1+2rz,2L^2,1];J=Mod[1+ 1+2ry,2L^2,1];*)
(* {{     1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0,0,0},{0,0,u[[ Mod[rz+1,L^2,1],2]] ,0},{0,0,0,0}}        ],{n,0,L-1} , {m,0,L-1}];*)
(*SS        ];*)
(**)


(* ::Input:: *)
(*(* writting the exp values in the form of the MF parameters:   *)*)
(*toMFparametersPure[U_,u_,L_,M_] :=Module[  { Nc=L^2,\[DoubleStruckCapitalU],\[DoubleStruckCapitalU]h,icc,SS=Array[Null,3]  },*)
(*\[DoubleStruckCapitalU]=TmatPure[L] . U;\[DoubleStruckCapitalU]h=\[DoubleStruckCapitalU][[;;,-Nc;;-1]];icc=Chop[I  \[DoubleStruckCapitalU]h . \[DoubleStruckCapitalU]h\[ConjugateTranspose] ];*)
(**)
(*SS[[3]]=ArrayFlatten[ # ,1]& @ Table[ *)
(*Module[{ rz,Io,J},rz=m+n L; Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1 +2rz,2L^2,1];*)
(* {{   -u[[rz+1,3]]  1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,-1 }}    ],{n,0,L-1} , {m,0,L-1}];*)
(**)
(*SS[[1]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,rx,mx,Io,J},mx=Mod[m+1,L]; rz=m+n L;rx= mx+n L ;Io=Mod[1+2rz,2L^2,1];J=Mod[1 + 1+2rx,2L^2,1];*)
(* {{ -u[[ rz+1 ,1]]  1/2 (icc[[J,Io]]-icc[[Io,J]] ) ,0,0,0},{0,-1 ,0,0},{0,0,0,0},{0,0,0,0}}     ] ,{n,0,L-1} , {m,0,L-1}];*)
(*SS[[2]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,ry,my,ny,Io,J}, ny=Mod[n+1,L];my=If[ny==0,Mod[m-M,L] , m];  ry= my+ny L;*)
(*rz=m+n L;Io=Mod[1+2rz,2L^2,1];J=Mod[1+ 1+2ry,2L^2,1];*)
(* {{   -u[[rz+1,2]]  1/2 (icc[[J,Io]]-icc[[Io,J]] ),0,0,0},{0,0,0,0},{0,0,-1 ,0},{0,0,0,0}}        ],{n,0,L-1} , {m,0,L-1}];*)
(*SS        ];*)
(**)


(* ::Subsubsection::Bold::Closed:: *)
(*gauge pure*)


(* ::Code::Bold:: *)
(**)


nx={1/2,Sqrt[3]/2};ny={-(1/2),Sqrt[3]/2};     \[Delta]x={-Sqrt[3]/2,-1/2}/Sqrt[3];\[Delta]y={Sqrt[3]/2,-1/2}/Sqrt[3];\[Delta]z={0,1}/Sqrt[3];
 
asites[m_,n_]:=m nx+n ny;
bsites[m_,n_]:=m nx+n ny-\[Delta]z;
gaugeConfigUPure[\[Sigma]_,L1_,L2_]:=Module[{bonds,linethick=.012 (Log[49]/Log[L1 L2]),lattice,size}, size = 0.01 (Log[49]/Log[L1 L2])^2 ;
lattice=Flatten[Table[{  {Black,PointSize[4 size],Point[asites[m,n]]},{Black,PointSize[4 size],Point[bsites[m,n]]}  , {Gray,PointSize[2  size ],Point[asites[m,n]]},{White,PointSize[2  size],Point[bsites[m,n]]}     },{m,0,L1-1},{n,0,L2-1} ], 2    ];
bonds=Flatten[Table[
{    
{ If[\[Sigma][[m+n L1+1,3]]>= 0, Cyan ,Darker@Blue],Thickness[linethick],Line[{asites[m,n],bsites[m,n]}]},                   (*z-bonds*)
{If[\[Sigma][[m+n L1+1,1]]>= 0, Magenta ,Darker@Red],Thickness[linethick],Line[{asites[m,n],bsites[m+1,n]}]},             (*x-bonds*)
{If[\[Sigma][[m+n L1+1,2]]>= 0, Darker@Yellow ,Darker@Green],Thickness[linethick],Line[{asites[m,n],bsites[m,n+1]}]}              (*y-bonds*)
},{m,0,L1-1},{n,0,L2-1} ], 2    ];
Graphics[{bonds,lattice},ImageSize->300]            ];  

uniformU[u_,L1_,L2_]:=Table[   {u,u,u}     ,L1 L2];
flipUoneZbond[u0_,L1_,L2_]:= Module[{R1,u=u0,m0=\[LeftCeiling]L1/2\[RightCeiling] -1,n0=\[LeftCeiling](L2+1)/2\[RightCeiling]-1 }, R1= { m0,n0};
  Module[{r,m,n}, m =Mod[R1[[1]],L1] ; n=Mod[R1[[2]],L2] ; r = m+n L1+1;    u[[r, 3 ]] =-u0[[r,3 ]];      ]   ;  u           ];
flipUZMaxSpaced[u0_,L1_,L2_]:= Module[{R1,d,u=u0,m0=\[LeftCeiling]L1/2\[RightCeiling] -1,n0=\[LeftCeiling](L2+1)/2\[RightCeiling]-1 }, d=\[LeftFloor](L1-1)/2\[RightFloor]; R1= { m0,n0};
Do[  Module[{r,m,n}, m =Mod[R1[[1]]+i,L1] ; n=Mod[R1[[2]]-i,L2] ; r = m+n L1+1;    u[[r, 3 ]] =-u0[[r,3 ]];      ]   , {i,0,d-1}];  u           ];


(* ::Text::Bold:: *)
(*		consider other gauge configurations: *)


(* ::Code::Bold:: *)
(**)


uniformU[u_,L_]:=Table[ {u,u,u},L^2];
gauge2vX[u0_,L_]:= Module[{d1,d2,u=u0,mS,nS,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ; r = m+n L+1;    u[[r,1]] =-u0[[r,1]];      ]   , {i,1,d1}];               u           ];
gauge2vY[u0_,L_]:= Module[{d1,d2,u=u0,mW,nW,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mW=r2-1;nW=L-r2-1;
Do[  Module[{r,m,n}, m =mW+i; n=nW; r = m+n L+1;    u[[r,2]] =-u0[[r,2]];      ]   , {i,1,d1}];                u           ];
gauge2vZ[u0_,L_]:= Module[{d1,d2,u=u0,mW,nW,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mW=r2-1;nW=L-r2-1;
Do[  Module[{r,m,n}, m =mW+i; n=nW-i+1; r = m+n L+1;    u[[r,3]] =-u0[[r,3]];      ]   , {i,1,d1}];                u           ];
gauge4v[u0_,L_]:= Module[{d1,d2,u=u0,mS,nS,mN,nN,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;       r = m+n L+1;    u[[r,1]] =-u0[[r,1]];      ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1; r = m+n L+1;    u[[r,1]] =-u0[[r,1]];      ]   , {i,1,d1}];               u           ];
positionVortex[v_,L_]:=Module[{d1,d2,r2,mS,nS,mW,nW,mE,nE,mN,nN},
		 r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2; 

		 mS=r2-1;           nS=r2-1;
		 mW=r2-1;           nW=L-r2-1;
		 mE=L-r2-1;     nE=r2-1;
		 mN=L-r2-1;     nN=L-r2-1;

		{{mS,nS},{mW,nW},{mE,nE},{mN,nN}}[[v]]
] ;


(* ::Text::Bold:: *)
(*v=1  South, v=2 West, v=3 East, v=4 North*)


(* ::Code::Bold:: *)
(**)


gauge2v[u0_,L_,\[Alpha]_]:= Module[{d1,d2,u=u0,m0,n0,r2 ,v},   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor];   Do[  Module[{r,m,n},
							m={r2-1,             r2-1+i,              r2-1+i   }[[\[Alpha]]]; 
							n={r2-1+i,       L-r2-1,             L-r2-i   }[[\[Alpha]]];
 r = m+n L+1;    u[[r,\[Alpha]]] =-u0[[r,\[Alpha]]];      ]   , {i,1,d1}];                u           ];


(* ::Subsubsection::Bold::Closed:: *)
(*correlation*)


(* ::Code::Bold:: *)
(*spin-spin correlation*)


ssC[\[Chi]_,\[Omega]_,r_,L_,M_]:=Module[ {m,n,mx,my,ny,R},n=\[LeftFloor](r-1)/L\[RightFloor];m=r-1-n L;mx=Mod[m+1,L];ny=Mod[n+1,L];my=If[ny==0,Mod[m-M,L] , m];
 R={ mx+n L +1,my+ny L +1,m+n L +1};
Table[ \[Omega][[1,r,\[Alpha]+1,0+1]]   \[Omega][[2,R[[\[Gamma]]],\[Beta]+1,0+1]]   - \[Chi][[\[Gamma],r,\[Alpha]+1,\[Beta]+1]] \[Chi][[\[Gamma],r,0+1,0+1]]  +  \[Chi][[\[Gamma],r,\[Alpha]+1,0+1]] \[Chi][[\[Gamma],r,0+1,\[Beta]+1]]
,
{\[Gamma],1,3}, {\[Alpha],1,3},{\[Beta],1,3}]
];

ssC2[\[Chi]_,\[Omega]_,m_,n_,L_]:=Module[ {mx,ny,R,r=Mod[m,L] +Mod[n,L] L +1},mx=Mod[m+1,L];ny=Mod[n+1,L]; R={ mx+n L +1,m+ny L +1,r};
Table[ \[Omega][[1,r,\[Alpha]+1,0+1]]   \[Omega][[2,R[[\[Gamma]]],\[Beta]+1,0+1]]   - \[Chi][[\[Gamma],r,\[Alpha]+1,\[Beta]+1]] \[Chi][[\[Gamma],r,0+1,0+1]]  +  \[Chi][[\[Gamma],r,\[Alpha]+1,0+1]] \[Chi][[\[Gamma],r,0+1,\[Beta]+1]],
{\[Gamma],1,3}, {\[Alpha],1,3},{\[Beta],1,3}]                                      ];



(* ::Subsubsection::Bold::Closed:: *)
(*electric field*)


(* ::Code::Bold:: *)
(**)


uniformK[K_,L1_,L2_]:=Table[   {K,K,K}     ,L1 L2];

addVortex[ Ko_,K0_, R_,L1_,L2_]:= Module[ {m=R[[1]],n=R[[2]],K=Ko},
K[[ m+n L1 +1,3    ]]=K0;
K[[ Mod[m+1,L1]+Mod[n-1,L2] L1 +1  ,2    ]]=K0;
K[[ Mod[m+1,L1]+n L1 +1,1    ]]=K0;
K[[ Mod[m+1,L1]+Mod[n+1,L2] L1 +1  ,3   ]]=K0;
K[[ m+Mod[n+1,L2] L1+1   ,2    ]]=K0;
K[[ Mod[m-1,L1]+Mod[n+1,L2] L1 +1  ,1   ]]=K0;
K];

add2Vortices[ Ko_,K1_,K2_, R1_,R2_,L1_,L2_] :=Module[ {K}, K = addVortex[ Ko,K1, R1,L1,L2]; addVortex[ K,K2, R2,L1,L2]];
add4Vortices[ Ko_,Kmod_, R1_,R2_, R3_,R4_,L_] :=Module[ {K}, K = addVortex[ Ko,Kmod, R1,L,L]; K = addVortex[ K,Kmod, R2,L,L]; K = addVortex[ K,Kmod, R3,L,L]; addVortex[ K,Kmod, R4,L,L]];

add2VorticesMaxSpaced[ Ko_,K1_,K2_,L1_,L2_] := Module[{R1,R2,d}, d=\[LeftFloor](L1-1)/2\[RightFloor]; R1= {\[LeftFloor]L1/2\[RightFloor] ,\[LeftFloor](L2+1)/2\[RightFloor] };R2= {\[LeftFloor]L1/2\[RightFloor]   +d,\[LeftFloor](L2+1)/2\[RightFloor] -d};
	add2Vortices[ Ko,K1,K2, R1,R2,L1,L2]];

add4VorticesMaxSpaced[ K0_,Kmod_,L_] :=Module[ {RS,RE,RW,RN,r2},   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor];
		RS={(r2-1),(r2-1)};
		RW={(r2-1),(L-r2-1)};
		RE={(L-r2-1),(r2-1)};
		RN={(L-r2-1),(L-r2-1)};
add4Vortices[ K0,Kmod, RS,RW, RE,RN,L] ];

\[Chi]Gx = {{0.5249,0,0,0},{0,-1,0,0},{0,0,0,0},{0,0,0,0}};\[Chi]Gy = {{0.5249,0,0,0},{0,0,0,0},{0,0,-1,0},{0,0,0,0}};\[Chi]Gz ={{0.5249,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,-1}};

\[Omega]GA = {{I,0.000001,0.000001,0.000001},{-0.000001,I,0.000001,-0.000001},{-0.000001,-0.000001,I,0.000001},{-0.000001,0.000001,-0.000001,I}}; \[Omega]GB = \[Omega]GA;



(* ::Subsubsection::Bold::Closed:: *)
(*more*)


(* ::Code::Bold:: *)
(**)


asites[m_,n_]:=m nx+n ny;
bsites[m_,n_]:=m nx+n ny-\[Delta]z;


(* ::Subsection:: *)
(*Momentum*)


(* ::Subsubsection::Closed:: *)
(*consts*)


mx= 2\[Pi] {1 , 1/Sqrt[3]}; my= 2\[Pi] {-1 , 1/Sqrt[3]};
(*toMomentum[n_,Nb_]:= If[FractionalPart[Sqrt@Nb] \[NotEqual] 0, Print["Nb=",Nb,", is not a perfect square."];Abort[],
				Module[ {L = IntegerPart[Sqrt@Nb],p1,p2},  p1 =  Mod[n-1,L]+1; p2 =   (n-p1)/L+1;   mxp1/L +myp2/L]      ];*)
toMomentum[n_,Nb_]:= Module[ {L = IntegerPart[Sqrt@Nb],p1,p2},  p1 =  Mod[n,L,1]; p2 =   (n-p1)/L+1;   mx p1/L +my p2/L];
toMomentumTable[L_]:=Table[Module[{p1,p2},p1=Mod[n,L,1];p2=(n-p1)/L+1;{(2(p1-p2)\[Pi])/L,(2 (p1+p2) \[Pi])/(Sqrt[3] L)}],{n,1,L^2}];
toMomentumTable[L_]:=Table[Module[{m,n},m=Mod[r-1,L];n=\[LeftFloor](r-1)/L\[RightFloor];{(2\[Pi](m-n))/L,(2\[Pi] (m+n))/(Sqrt[3] L)}],{r,1,L^2}];
toMomentumInverse[k_,Nb_]:= If[FractionalPart[Sqrt@Nb] != 0, Print["Nb=",Nb,", is not a perfect square."];Abort[],    
Module[ {M = IntegerPart[Sqrt@Nb],p1,p2},   p1=k . nx M/(2\[Pi]); p2=k . ny M/(2\[Pi]);  Round[p1 + (p2-1)M ]  ]];
kD2=1/3 mx+2/3 my;
kD1=2/3 mx+1/3 my;
\[CapitalGamma]point={0,0};
Kpoint=2/3 mx+1/3 my;
KPpoint=1/3 mx+2/3 my;
Mpoint=1/2 (Kpoint+KPpoint);
MPpoint=1/2 (Kpoint-mx+KPpoint);
MPPpoint=1/2 (Kpoint+KPpoint-my);


toLine\[CapitalGamma]KK\[CapitalGamma]0[M_]:=Module[ {ki=my,kf=mx},  Table[ ki+t(kf-ki),{t,0,1,1/M}]    ] ;
toLine\[CapitalGamma]KM\[CapitalGamma][M_]:=Module[ {ks,Ne,\[CapitalDelta]t},(* 1: \[CapitalGamma]=0, 2: K=(1/3,2/3), 3: K'=(2/3,3/3), 4: \[CapitalGamma]=0,  *)
ks={ \[CapitalGamma]point,Kpoint,MPpoint,\[CapitalGamma]point }; 
Ne=Length[ks]-1; (* number of edges making the path in the BZ*) 
Flatten[ Table[ Table[ ks[[i]]+t(  ks[[i+1]]-ks[[i]]  ),{t,0,1-Ne/M (1-KroneckerDelta[i,Ne]),Ne/M}],{i,1,Ne}],1]       ] ;
toLine\[CapitalGamma]KM\[CapitalGamma][M_]:=Module[ {ks,Ne,\[CapitalDelta]t},(* 1: \[CapitalGamma]=0, 2: K=(1/3,2/3), 3: K'=(2/3,3/3), 4: \[CapitalGamma]=0,  *)
ks={ \[CapitalGamma]point,Kpoint,MPpoint,\[CapitalGamma]point }; 
Ne=Length[ks]-1; (* number of edges making the path in the BZ*) 
Flatten[ Table[ Table[ ks[[i]]+t(  ks[[i+1]]-ks[[i]]  ),{t,Ne/M (1-KroneckerDelta[i,1]),1,Ne/M}],{i,1,Ne}],1]      ] ;
toLine\[CapitalGamma]KM\[CapitalGamma]KMM[M_]:=Module[ {ks,Ne,\[CapitalDelta]t},(* 1: \[CapitalGamma]=0, 2: K=(1/3,2/3), 3: K'=(2/3,3/3), 4: \[CapitalGamma]=0,  *)
ks={ \[CapitalGamma]point,KPpoint,MPpoint,\[CapitalGamma]point,Kpoint,Mpoint,MPPpoint }; 
Ne=Length[ks]-1; (* number of edges making the path in the BZ*) 
Flatten[ Table[ Table[ ks[[i]]+t(  ks[[i+1]]-ks[[i]]  ),{t,Ne/M (1-KroneckerDelta[i,1]),1,Ne/M}],{i,1,Ne}],1]
    ] ;
    


(* ::Subsubsection::Closed:: *)
(*ham*)


(*Tmom[k_]:=Module[{kx=k.mx,ky=k.my},{{1,0,0,0,1,0,0,0},{0,1,0,0,0,1,0,0},{0,0,1,0,0,0,1,0},{0,0,0,1,0,0,0,1},{\[ImaginaryI],0,0,0,-\[ImaginaryI],0,0,0},{0,\[ImaginaryI] \[ExponentialE]^(\[ImaginaryI] kx),0,0,0,-\[ImaginaryI] \[ExponentialE]^(\[ImaginaryI] kx),0,0},{0,0,\[ImaginaryI] \[ExponentialE]^(\[ImaginaryI] ky),0,0,0,-\[ImaginaryI] \[ExponentialE]^(\[ImaginaryI] ky),0},{0,0,0,\[ImaginaryI],0,0,0,-\[ImaginaryI]}}];*)
Tkmom=KroneckerProduct[ {{1,1},{I,-I}}, SparseArray[{Band[{1,1}]-> {1,1,1,1}},{4,4}] ]; 
AMFmom[Jm_,U_]:=Module[{AMF,N=Nmat},  AMF=Table[    Sum[   2 Jm[[\[Gamma]]][[\[Alpha],\[Beta]]] N[[\[Alpha]]] . U[[\[Gamma]]] . N[[\[Beta]]] ,{\[Alpha],1,3},{\[Beta],1,3}]    ,{\[Gamma],1,3}];   KroneckerProduct[ {{0,1},{0,0}},#]&/@Re@AMF];
BMFmom[Jm_,h_,\[Lambda]_,V_]:=Module[{Ba,Bb,M=Mmat,N=Nmat,G=Gmat},
Ba=Sum[    Sum[    Jm[[\[Gamma]]][[\[Alpha],\[Beta]]] N[[\[Alpha]]] traceN[ V[[2]] ][[\[Beta]]] ,{\[Alpha],1,3},{\[Beta],1,3}] -2 h[[\[Gamma]]] N[[\[Gamma]]] + \[Lambda][[1]][[\[Gamma]]] G[[\[Gamma]]]        ,{\[Gamma],1,3}];
Bb=Sum[    Sum[    Jm[[\[Gamma]]][[\[Alpha],\[Beta]]] N[[\[Alpha]]] traceN[ V[[1]] ][[\[Beta]]] ,{\[Alpha],1,3},{\[Beta],1,3}] -2 h[[\[Gamma]]] N[[\[Gamma]]] + \[Lambda][[2]][[\[Gamma]]] G[[\[Gamma]]]         ,{\[Gamma],1,3}];
{KroneckerProduct[ {{1,0},{0,0}},Re@Ba],KroneckerProduct[ {{0,0},{0,1}},Re@Bb] }];

cMFmom[Jm_,U_,V_]:=Module[{N=Nmat}, Sum[ 1/8  Jm[[\[Gamma]]][[\[Alpha],\[Beta]]] (traceN[ V[[1]] ][[\[Alpha]]]  traceN[ V[[2]] ][[\[Beta]]] +2  Tr[ U[[\[Gamma]]]\[Transpose] . N[[\[Alpha]]] . U[[\[Gamma]]] . N[[\[Beta]]] ]   ),{\[Alpha],1,3},{\[Beta],1,3},{\[Gamma],1,3}]          ];
enMFmom[Jm_,U_,V_,h_,\[Eta]_:0]:=Module[{M=Mmat,N=Nmat,G=Gmat,\[Lambda]},\[Lambda]=\[Eta]{h,h}+\[Lambda]effmom[Jm,h,V];  cMFmom[Jm,U,V]+1/4  Sum[-2h[[\[Gamma]]] traceN[ V[[\[Sigma]]] ][[\[Gamma]]] + \[Lambda][[\[Sigma],\[Gamma]]] traceG[ V[[\[Sigma]]] ][[\[Gamma]]]
,{\[Sigma],1,2},{\[Gamma],1,3}]        ];
enSUMmom[Jm_,U_,V_,h_,L_:30,\[Eta]_:0]:=Module[{\[Lambda],mT=toMomentumTable[L]},-cMFmom[Jm,U,V]+1/(2L^2) Sum[Total@Select[Eigenvalues@N@HmfMomentum[Jm,h,U,V,mT[[i]],\[Eta]],#<=0& ], {i,1,L^2}] ];

HmfMomentum[Jmatrice_,h_,U_,V_, k_,\[Eta]_:0] := Module[{Ax,Ay,Az,BA,BB, hx,hy,hz, kx,ky,ha,HA,HB,\[Lambda]}, kx=k . nx;ky=k . ny;
\[Lambda]=\[Eta]{h,h}+ \[Lambda]effmom[Jmatrice,h,V];
{Ax,Ay,Az}=I AMFmom[Jmatrice,U];  ha= Exp[ I kx]  Ax + Exp[ I ky] Ay + Az;   
 (*ha=Exp[-I k.\[Delta]x] Ax + Exp[-I k.\[Delta]y] Ay +Exp[-I k.\[Delta]z] Az;*)   HA=N@( ha+ConjugateTranspose@ha );
{BA,BB}= I BMFmom[Jmatrice,h,\[Lambda],V];   HB=BA+BB ; (*  HB=N[   ( HB+ConjugateTranspose@HB )/2 ];*)
1/2 (HA+HB)   ];
UmatK[H_]:= Module[ {R=Eigensystem@N[H]},ReverseSort[ R\[Transpose] ]\[Transpose][[ 2 ]]\[Transpose] ];
\[Lambda]effmom[Jmat_,h_,V_]:=Module[{M},M=Table[1/2 traceN[V[[\[Sigma]]] ],{\[Sigma],1,2} ];  Table[ -h + 1/2 Sum[ Jmat[[\[Gamma]]] . M[[  Mod[\[Sigma]+1,2,1]  ]]  ,{\[Gamma],1,3}],    {\[Sigma],1,2}]     ];
HmfMomentumVec[Jmatrice_,h_,U_,V_, kTable_,\[Eta]_:0] :=Table[HmfMomentum[Jmatrice,h,U,V, kTable[[l]],\[Eta]],{l,1,Length@kTable}];


UmatVec[Jmatrice_,h_,U_,V_, kTable_,Tk_,\[Eta]_:0] :=UmatK/@Table[ Tk\[ConjugateTranspose] . HmfMomentum[Jmatrice,h,U,V, kTable[[l]],\[Eta]] . Tk,{l,1,Length@kTable}];


(* ::Subsection::Bold:: *)
(*MF model definitions*)


(* ::Subsubsection::Bold::Closed:: *)
(*Saving and Loading data*)


cyclicPermutation[A_,s_:1]:= If[ Length[A]==3 \[And] Length[A[[1]]]==3, Table[ A[[Mod[i-s,3,1],Mod[j-s,3,1]]],{i,1,3},{j,1,3}]  ];
fromJmat[Jm_]:= Module[{J,K,\[CapitalGamma],\[CapitalGamma]p,DM,jm}, jm=(  cyclicPermutation[Jm[[1]] ,2] + cyclicPermutation[Jm[[2]] ,1] + Jm[[3]] )/3; 
\[CapitalGamma]p=(jm[[1,3]]+jm[[2,3]]+jm[[3,1]]+jm[[3,2]])/4; 
K=( jm[[3,3]]-(jm[[2,2]] +jm[[1,1]] )/2); J=(jm[[3,3]]-K);DM=(jm[[1,2]]-jm[[2,1]])/2;\[CapitalGamma]=(jm[[1,2]]+jm[[2,1]])/2;
N@{J,K,\[CapitalGamma],\[CapitalGamma]p,DM}    ];
createDir[path_] :=
Module[ {l=Length@FileNames[path]},		
	If[ l==0,  CreateDirectory@File@FileNameJoin[{path,"data" }]; CreateDirectory@File@FileNameJoin[{path ,"graph" }];  , Null]  ];
loadData[pathData_]:=
Module[{f,data},  f = OpenRead[pathData];
If[f==$Failed, Print["Failed to OpenRead file at: "]; Print[ pathData ]; Abort[]  ];
	data=ReadList[f]; Close[f]; data[[-1]]   ];
toPath800 [parameters0_,L_,acuracy_,gauge_,NbName_]:= Module[{h ,hS,parameters=parameters0, filename,r,\[Phi],\[Theta]},h=parameters[[2]] ;hS=parameters[[4]]; {r,\[Theta],\[Phi]}=hS;   
filename="t=X1_eV=X2_JKG=X3_JKGmod=X4";
If[gauge=="free", parameters[[3]]=parameters0[[1]]; parameters[[6]]=0; filename="free_X2_JKG=X3";];
If[gauge=="free_lambda",parameters[[3]]=parameters0[[1]];filename="lambda_X2_JKG=X3"; ];
FileNameJoin[{       Directory[]   , "Files" ,NbName,gauge,StringReplace[filename,{"X1"->  ToString[parameters[[5]]],"X2"->  ToString[parameters[[6]]],
"X3"->  ToString[NumberForm[#,{4,4}]&/@fromJmat@parameters[[1]]  ],
"X4"->  ToString[NumberForm[#,{4,4}]&/@fromJmat@parameters[[3]]]   }] , "data"  , 
StringReplace["h=(M,N,T)_L=Y_A=Z.txt",{"Y"-> ToString[L], "Z"-> ToString[acuracy],"M"->  ToString[r,InputForm] ,"N"->ToString@\[Theta],"T"->ToString@\[Phi]}   ] 
   }]];
dataToFile800[ parameters0_,L_,acuracy_,data_,gauge_,NbName_] :=
Module[ {path,f,filename,parameters=parameters0}, 
filename="t=X1_eV=X2_JKG=X3_JKGmod=X4";
If[gauge=="free",parameters[[3]]=parameters0[[1]];parameters[[6]]=0;filename="free_X2_JKG=X3";];
If[gauge=="free_lambda",parameters[[3]]=parameters0[[1]];filename="lambda_X2_JKG=X3"; ];
		createDir@FileNameJoin[{ Directory[] ,"Files",  NbName,gauge, StringReplace[filename,
{"X1"->ToString[parameters[[5]]],"X2"->ToString[parameters[[6]]],"X3"->ToString[NumberForm[#,{4,4}]&/@fromJmat@parameters[[1]]  ],"X4"->ToString[NumberForm[#,{4,4}]&/@fromJmat@parameters[[3]]]   }]     }] ;
		path = toPath800[parameters,L,acuracy,gauge,NbName];		
		f = OpenAppend[path];
		 Write[f, data];
		 Close[f];                ];
loadDataTry[pathData_]:=
Module[ {f,data,path=FindFile[pathData]},
If[path==$Failed,(*Print["New entry at:",pathData];*)Return[$Failed]];
		        f = OpenRead[pathData];
(*If[f==$Failed, Print["Failed to OpenRead file at: ", pathData ]; Abort[] ];*)
		        data=ReadList[f];
		        Close[f];			data[[-1]]
];


(* ::Subsubsection::Bold::Closed:: *)
(*Hamiltonian  matrices*)


Jmat[J_,K_,G_,Gp_,D_]:={
{{J[[1]]+K[[1]],Gp[[1]],Gp[[1]]},{Gp[[1]],J[[1]],G[[1]]+D[[1]]},{Gp[[1]],G[[1]]-D[[1]],J[[1]]}},
{{J[[2]],Gp[[2]],G[[2]]-D[[2]]},{Gp[[2]],J[[2]]+K[[2]],Gp[[2]]},{G[[2]]+D[[2]],Gp[[2]],J[[2]]}},
{{J[[3]],G[[3]]+D[[3]],Gp[[3]]},{G[[3]]-D[[3]],J[[3]],Gp[[3]]},{Gp[[3]],Gp[[3]],J[[3]]+K[[3]]}}   
 };
Hmf[Jmat_,h_,U_,V_,L_,\[Eta]_:1] := Module[{Nc=L^2,H,Hmat,\[Lambda]}, 
\[Lambda]=\[Eta] \[Lambda]eff[Jmat,h,V,L];
Hmat=Haux[Jmat,h,U,V,\[Lambda],L];
H= I Sum[Module[{rx,ry,rz,mx,ny}, mx=Mod[m+1,L];ny=Mod[n+1,L];
 rx=mx+n L+1;ry=m+ny L+1;rz=m+n L+1;
  KroneckerProduct[ one[rz,rx,Nc],Hmat[[rz,1]]      ]+ 
KroneckerProduct[ one[rz,ry,Nc],Hmat[[rz,2]]    ]+ 
KroneckerProduct[ one[rz,rz,Nc],Hmat[[rz,3]]    ]
] ,{m,0,L-1},{n,0,L-1}];  1/2 (H+H\[ConjugateTranspose])];
Haux[Jmat_,h_,U_,V_,\[Lambda]_,L_]:=Module[{M=Mmat,Nm=Nmat,G=Gmat},Table[Module[{ Ax,Ay,Az,BA,BB},
{Ax,Ay,Az}=I KroneckerProduct[ {{0,1},{0,0}},#]&/@Re@Table[    Sum[   2  Jmat[[r,\[Gamma]]][[\[Alpha],\[Beta]]] Nm[[\[Alpha]]] . U[[r,\[Gamma]]] . Nm[[\[Beta]]] ,{\[Alpha],1,3},{\[Beta],1,3}], {\[Gamma],1,3}]   ;
{BA,BB}        =I{
KroneckerProduct[ {{1,0},{0,0}},Re@Sum[Sum[ Jmat[[r,\[Gamma]]][[\[Alpha],\[Beta]]] Nm[[\[Alpha]]] Tr[Transpose[V[[r,2]]] . Nm[[\[Beta]]] ],{\[Alpha],1,3},{\[Beta],1,3}] -2 h[[\[Gamma]]] Nm[[\[Gamma]]]+\[Lambda][[r,1]][[\[Gamma]]] G[[\[Gamma]]],{\[Gamma],1,3}]],
KroneckerProduct[ {{0,0},{0,1}},Re@Sum[Sum[ Jmat[[r,\[Gamma]]][[\[Alpha],\[Beta]]] Nm[[\[Alpha]]] Tr[Transpose[V[[r,1]]] . Nm[[\[Beta]]] ],{\[Alpha],1,3},{\[Beta],1,3}] -2 h[[\[Gamma]]] Nm[[\[Gamma]]]+\[Lambda][[r,2]][[\[Gamma]]] G[[\[Gamma]]],{\[Gamma],1,3}]]   };
{Ax,Ay,Az+BA+BB} ],{r,1,L^2}]                  ];
UmatK[H_]:= Module[ {R=Eigensystem@N[H]},ReverseSort[R\[Transpose]]\[Transpose][[2]]\[Transpose] ];
\[Lambda]eff[Jmat_,h_,V_,L_]:=Module[{M},  M=Table[1/2 traceN[V[[r,\[Sigma]]] ],{r,1,L^2},{\[Sigma],1,2} ] ; Table[ h-1/2 Sum[ 
Jmat[[r,\[Gamma]]] . M[[   Mod[r+{ Mod[r,L]-Mod[r-1(2\[Sigma]-3),L] ,(Mod[\[LeftFloor](r-1)/L\[RightFloor]+1(2\[Sigma]-3),L]-Mod[\[LeftFloor](r-1)/L\[RightFloor],L])  L,0 }[[\[Gamma]]]  ,L^2,1],Mod[\[Sigma]+1,2,1]]]  ,{\[Gamma],1,3}],{r,1,L^2},{\[Sigma],1,2}]    ];


(* ::Subsubsection::Bold::Closed:: *)
(*Mean field parameters *)


Tx = KroneckerProduct[ {{1,1},{I,-I}}, SparseArray[{Band[{1,1}]-> {0,1,0,0}},{4,4}]  ]; 
Ty = KroneckerProduct[ {{1,1},{I,-I}}, SparseArray[{Band[{1,1}]-> {0,0,1,0}},{4,4}]  ];
Tz = KroneckerProduct[ {{1,1},{I,-I}}, SparseArray[{Band[{1,1}]-> {1,0,0,1}},{4,4}]  ]; 
Tmat[L1_,L2_] := Module[{Nc=L1 L2,T}, T=Sum[Module[{rx,ry,rz,mx,ny},
 mx=Mod[m+1,L1,1];ny=Mod[n+1,L2,1]; rx=mx+(n-1)L1;ry=m+(ny-1)L1;rz=m+(n-1)L1;
 KroneckerProduct[ one[rz,rx,Nc],Tx]+  KroneckerProduct[ one[rz,ry,Nc],Ty]+ KroneckerProduct[ one[rz,rz,Nc],Tz ]] ,{m,1,L1},{n,1,L2}];T];

bilinears[U_,L_] :=Module[  { Nc=L^2,TU,TUh,u,\[Chi]=Array[Null,3],\[Omega]=Array[Null,2],\[Xi]=Array[Null,6]},
TU=Tmat[L,L] . U;TUh=TU[[;;,-4Nc;;-1]];u=I  TUh . TUh\[ConjugateTranspose] ;
\[Chi][[3]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,Io,J},rz=m+n L; Io=Mod[\[Alpha]+8rz,8L^2,1];J=Mod[\[Beta] + 4 +8rz,8L^2,1];
1/2 ( u[[J,Io]] -u[[Io,J]] )   ],{n,0,L-1}, {m,0,L-1},  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Chi][[1]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,rx,mx,Io,J},mx=Mod[m+1,L ]; rz=m+n L;rx=mx+n L;Io=Mod[\[Alpha]+8rz,8L^2,1];J=Mod[\[Beta] + 4+8rx,8L^2,1];
1/2 ( u[[J,Io]] -u[[Io,J]] )    ]   ,{n,0,L-1}, {m,0,L-1},  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Chi][[2]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,ry,ny,Io,J}, ny=Mod[n+1,L];  ry=m+ny L;rz=m+n L;Io=Mod[\[Alpha]+8rz,8L^2,1];J=Mod[\[Beta] + 4+8ry,8L^2,1];
1/2 ( u[[J,Io]] -u[[Io,J]])  ]    , {n,0,L-1}, {m,0,L-1},{\[Alpha],1,4},{\[Beta],1,4} ];
\[Omega][[1]]= ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,Io,J},rz=m+n L; Io=Mod[\[Alpha]+8rz,8L^2,1];J=Mod[\[Beta] +8rz,8L^2,1];
I KroneckerDelta[Io,J]  + 1/2 ( u[[J,Io]] -u[[Io,J]]) (*+ (1/2)(u[[J,Io]] -u[[Io,J]])  u[[J,Io]]  *)   ], {n,0,L-1}, {m,0,L-1},{\[Alpha],1,4},{\[Beta],1,4} ];
\[Omega][[2]]= ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,Io,J},rz=m+n L; Io=Mod[\[Alpha]+4+8rz,8L^2,1];J=Mod[\[Beta] +4+8rz,8L^2,1];
I KroneckerDelta[Io,J] +1/2 ( u[[J,Io]] -u[[Io,J]]) (*+ (1/2)(u[[J,Io]] -u[[Io,J]])  u[[J,Io]]  *)  ], {n,0,L-1}, {m,0,L-1},  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Xi][[1]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,r,Io,J}, rz=m+n L;r=m+Mod[n+1,L ] L;Io=Mod[\[Alpha]+8rz,8L^2,1];J=Mod[\[Beta] +8r,8L^2,1];
1/2 ( u[[J,Io]] -u[[Io,J]])         ]    , {n,0,L-1}, {m,0,L-1},{\[Alpha],1,4},{\[Beta],1,4} ];
\[Xi][[2]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,r,Io,J}, rz=m+n L;r=Mod[m-1,L ]+n L;Io=Mod[\[Alpha]+8rz,8L^2,1];J=Mod[\[Beta] +8r,8L^2,1];
1/2 ( u[[J,Io]] -u[[Io,J]])           ]   ,{n,0,L-1}, {m,0,L-1},  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Xi][[3]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,r,Io,J}, rz=m+n L;  r=Mod[m+1,L ]+Mod[n-1,L] L; Io=Mod[\[Alpha]+8rz,8L^2,1];J=Mod[\[Beta] +8r,8L^2,1];
1/2 ( u[[J,Io]] -u[[Io,J]])            ]    , {n,0,L-1}, {m,0,L-1},{\[Alpha],1,4},{\[Beta],1,4} ];
\[Xi][[4]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,r,Io,J},  rz=m+n L;r=m+Mod[n-1,L ] L;Io=Mod[\[Alpha]+ 4+8rz,8L^2,1];J=Mod[\[Beta] + 4+8r,8L^2,1];
1/2 ( u[[J,Io]] -u[[Io,J]])            ]    , {n,0,L-1}, {m,0,L-1},{\[Alpha],1,4},{\[Beta],1,4} ];
\[Xi][[5]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,r,Io,J},  rz=m+n L;r=Mod[m+1,L ]+n L;Io=Mod[\[Alpha]+ 4+8rz,8L^2,1];J=Mod[\[Beta] + 4+8r,8L^2,1];
1/2 ( u[[J,Io]] -u[[Io,J]])            ]   ,{n,0,L-1}, {m,0,L-1},  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Xi][[6]]=ArrayFlatten[ # ,1]& @ Table[ Module[{ rz,r,mx,ny,Io,J},
r=Mod[m-1,L ]+Mod[n+1,L] L;rz=m+n L; Io=Mod[\[Alpha]+4+8rz,8L^2,1];J=Mod[\[Beta]+4 +8r,8L^2,1];
1/2 ( u[[J,Io]] -u[[Io,J]])           ]    , {n,0,L-1}, {m,0,L-1},{\[Alpha],1,4},{\[Beta],1,4} ];
Chop[{\[Chi][[1]],\[Chi][[2]],\[Chi][[3]],\[Omega][[1]],\[Omega][[2]],{\[Xi][[1]],\[Xi][[2]],\[Xi][[3]]},{\[Xi][[4]],\[Xi][[5]],\[Xi][[6]]} ,u },10^-16]
];
\[Chi]gauge4v[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;       r = m+n L+1;    \[Chi][[1,r]][[2,2]] =-\[Chi]0[[1,r]][[2,2]];      ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1; r = m+n L+1;    \[Chi][[1,r]][[2,2]] =-\[Chi]0[[1,r]][[2,2]];       ]   , {i,1,d1}];       \[Chi] ];
\[Chi]fixgauge[\[Chi]0_,u0_,L_]:= Module[{\[Chi]=\[Chi]0 }, 
Do[  \[Chi][[\[Alpha],r]][[\[Alpha]+1,\[Alpha]+1]] = u0[[r,\[Alpha]]]Abs@\[Chi]0[[\[Alpha],r]][[\[Alpha]+1,\[Alpha]+1]];      ,{\[Alpha],1,3} , {r,1,L^2}];   
   \[Chi] ];
icc[U_,L_,T_]:=Module[  { Nc=L^2,TUh,icc}, TUh=Chop[(T . U)[[;;,-4Nc;;-1]],10^-12];icc= I TUh . TUh\[ConjugateTranspose] ; (icc-ConjugateTranspose[icc] )/2 ];



(* ::Subsubsection::Bold::Closed:: *)
(*Energy MF*)


constantMF[Jmat_,U_,V_,L_]:= Module[{Nc=L^2,C},
C=Sum[Module[{r,rA,rB,mx,ny,N,G}, 	
mx=Mod[m+1,L];ny=Mod[n+1,L];
r=m+n L+1;   
rA={mx+n L+1,m+ny L+1,m+n L+1};
mx=Mod[m-1,L];ny=Mod[n-1,L];
rB={mx+n L+1,m+ny L+1,m+n L+1};
 1/8 Sum[ traceN[V[[rA[[3]],1]]] . Jmat[[rA[[3]],\[Gamma]]] . traceN[V[[rA[[\[Gamma]]],2]]]+ 
traceN[V[[rB[[3]],2]]] . Jmat[[rB[[\[Gamma]]],\[Gamma]]] . traceN[V[[rB[[\[Gamma]]],1]]],{\[Gamma],1,3}]
+1/4 Sum[    2 Tr[Transpose[Jmat[[r,\[Gamma]]]] . traceUNUN[ U[[r,\[Gamma]]] ]]  ,{\[Gamma],1,3}] 
    ],{m,0,L-1},{n,0,L-1} ];   Re@(C/( 2Nc ))   ];
EnMF[Jmat_,h_,U_,V_,L_,\[Eta]_:1]:= Module[{Nc=L^2,E,\[Lambda]},
If[\[Eta]==1,\[Lambda]=\[Lambda]eff[Jmat,h,V,L],\[Lambda]=\[Eta] \[Lambda]eff[Jmat,h,V,L]  ];

E=Sum[Module[{r,rA,rB,mx,ny,N,G}, 	
mx=Mod[m+1,L];ny=Mod[n+1,L];
r=m+n L+1;   
rA={mx+n L+1,m+ny L+1,m+n L+1};
mx=Mod[m-1,L];ny=Mod[n-1,L];
rB={mx+n L+1,m+ny L+1,m+n L+1};

1/4 ( \[Lambda][[r,1]] . traceN[V[[r,1]]]+\[Lambda][[r,2]] . traceN[V[[r,2]]]  ) 
-1/2 h . (   traceN[V[[r,1]]]+traceN[V[[r,2]]]  )         
+ 1/8 Sum[ traceN[V[[rA[[3]],1]]] . Jmat[[rA[[3]],\[Gamma]]] . traceN[V[[rA[[\[Gamma]]],2]]]+ 
traceN[V[[rB[[3]],2]]] . Jmat[[rB[[\[Gamma]]],\[Gamma]]] . traceN[V[[rB[[\[Gamma]]],1]]],{\[Gamma],1,3}]
+1/4 Sum[    2 Tr[Transpose[Jmat[[r,\[Gamma]]]] . traceUNUN[ U[[r,\[Gamma]]] ]]  ,{\[Gamma],1,3}] 
    ],{m,0,L-1},{n,0,L-1} ]; 

Re@(E/( 2Nc ))   ];
eigenvaluesEmf[Jmat_,h_,U_,V_,L_,\[Lambda]_]:=-constantMF[Jmat,U,V,L]+1/(2L^2) Total@Select[Quiet@Eigenvalues@N@Hmf[Jmat,h,U,V,L,\[Lambda]] ,#<=0& ];


(* ::Subsubsection::Closed:: *)
(*Wilson loop*)


WilsonLoop[U_,V_,L_]:=Module[{cyc}, cyc={{4,3},{2,4},{3,2}};
Table[Module[{m,n,R,Wilson={0,0,0,0}},  n=\[LeftFloor]p/L\[RightFloor];m=p-n L; 
R={    (*  r,\[Sigma],\[Gamma] -> ext *)
{m+n L +1,0,3  },
{ Mod[m+1,L]+Mod[n,L] L +1  ,1,2 },
{Mod[m+1,L]+n L +1,0,1    },
{Mod[m+1,L]+Mod[n+1,L] L +1  ,1,3  },
{ m+Mod[n+1,L] L+1   ,0,2   },
{ Mod[m,L]+Mod[n+1,L] L +1  ,1,1 }
};
Wilson[[1]]= Product[   Extract[    V[[ R[[i]][[1]], R[[i,2]]+1  ]],cyc[[  R[[i,3]]  ]]    ]  ,{i,1,6}];
Wilson[[2]]= Product[   Module[{\[Gamma]=cyc[[R[[i,3]]]] [[g]]         },           U[[ R[[i]][[1]] , \[Gamma]-1  ]][[\[Gamma],\[Gamma]]]
]  ,{i,{1,3,5}},{g,1,2}];
(Wilson[[1]]+Wilson[[2]])
],{p,0,L^2-1}]  ];
setP[x_,p_]:=N[ Round[10^p x]/10^p];
setP1[x_]:=N[ Round[10 x]/10];
applyHueP[x_,p_]:=Hue@@Table[ setP[x[[n]],p[[n]] ],{n,1,4}]
WilsonLoopGraphicsZoom[U_,V_,L_]:=Module[{cyc,linethick=.01 (Log[7]/Log[L]) },(* cyc={{3,4},{4,2},{2,3}};*)cyc={{4,3},{2,4},{3,2}};
Flatten[#,1]&@Table[Module[{R,Wilson={0,0,0,0}},  
R={    (*  r,\[Sigma],\[Gamma] -> ext *)
{m+n L +1,0,3  },
{ Mod[m+1,L]+n L +1  ,1,2 },
{Mod[m+1,L]+n L +1,0,1    },
{Mod[m+1,L]+Mod[n+1,L] L +1  ,1,3  },
{ m+Mod[n+1,L] L+1   ,0,2   },
{ m+Mod[n+1,L] L +1  ,1,1 }
};
Wilson[[1]]= Product[   Extract[    V[[ R[[i]][[1]], R[[i,2]]+1  ]],cyc[[  R[[i,3]]  ]]    ]  ,{i,1,6}];
Wilson[[2]]= Product[   Module[{\[Gamma]=cyc[[R[[i,3]]]] [[g]]         },           U[[ R[[i]][[1]] , \[Gamma]-1  ]][[\[Gamma],\[Gamma]]]
]  ,{i,{1,3,5}},{g,1,2}];
{applyHueP[#,{3,1,3,1}] & @ Blend[{ Hue[0.97,0.3,0.96], Hue[0.5,0,0.5,0.7] ,Hue[0.65,1,1,0.4]},(1/2+1/2 (Wilson[[2]] ) )  ]   ,EdgeForm[{Thickness[5linethick],Hue[0.5,0.0,1,0.9]}  ], 
Polygon[{\[Sigma]sites[m+1,n,0],   \[Sigma]sites[m+1,n+1,1],\[Sigma]sites[m,n+1,0],\[Sigma]sites[m,n+1,1] ,\[Sigma]sites[m,n,0],\[Sigma]sites[m+1,n,1]}  ]
}  ],{m,0,L/2-1},{n,L/2,L-1}]  ];
MFgaugeConfigZoom[U0_,V_,L_,title_:"",comment_:{"",""}]:=Module[{bonds,linethick=.012 (Log[7]/Log[L]),lattice,size,U,rt=0.8,max,plaquettes}, size = 0.01 (Log[7]/Log[L])^2 ; max=Max@Abs@U0; U=U0/max;
lattice=Flatten[Table[{  {Black,PointSize[4 size],Point[asites[m,n]]},{Black,PointSize[4 size],Point[bsites[m,n]]}  , {Gray,PointSize[2  size ],Point[asites[m,n]]},{White,PointSize[2  size],Point[bsites[m,n]]}     },{m,0,L/2-1},{n,L/2,L-1} ], 2    ];
bonds=Flatten[Table[
{{ Blend[{ Cyan, Gray ,Darker@Magenta},1/2 (1-root[#,rt]&@U[[m+n L+1,3]][[4,4]])  ],Thickness[(0.5+0.5root[#,rt]&@Abs@U[[m+n L+1,3]][[4,4]])linethick],Line[{asites[m,n],bsites[m,n]}]},                   (*z-bonds*)
{ Blend[{ Cyan, Gray ,Darker@Magenta},1/2 (1-root[#,rt]&@U[[m+n L+1,1]][[2,2]])  ],Thickness[(0.5+0.5root[#,rt]&@Abs@\[Chi][[1,m+n L+1]][[2,2]])linethick],Line[{asites[m,n],bsites[m+1,n]}]},             (*x-bonds*)
{ Blend[{ Cyan, Gray ,Darker@Magenta},1/2 (1-root[#,rt]&@U[[m+n L+1,2]][[3,3]])  ],Thickness[(0.5+0.5root[#,rt]&@Abs@U[[m+n L+1,2]][[3,3]])linethick],Line[{asites[m,n],bsites[m,n+1]}]}              (*y-bonds*)
},{m,0,L/2-1},{n,L/2,L-1} ], 2    ];
plaquettes=WilsonLoopGraphicsZoom[U,V,L];
Graphics[{plaquettes,bonds,lattice},ImageSize->500,ContentSelectable->True,Epilog->{Inset[ MaTeX[title,Magnification->1.8] ,Scaled[{0,0.95}],{Left,Center}],Inset[ MaTeX[comment[[1]],Magnification->1.5] ,Scaled[{1.0,0.95}],{Right,Center}]
,Inset[ MaTeX[comment[[2]],Magnification->1.35] ,Scaled[{1.0,0.88}],{Right,Center}]}
]  ];



(* ::Subsubsection::Bold::Closed:: *)
(*Couplings  *)


Jc[eV0_,JH_,U_,t_]:=1/54 ((2 (t[[1]]-t[[3]])^2)/(eV0-JH+U)-(2 (t[[1]]-t[[3]])^2)/(eV0+JH-U)-(6 t[[1]] (t[[1]]+2 t[[3]]))/(eV0+3 JH-U)+(6 t[[1]] (t[[1]]+2 t[[3]]))/(eV0-3 JH+U)+(2 t[[1]]+t[[3]])^2/(-eV0+2 JH+U)+(2 t[[1]]+t[[3]])^2/(eV0+2 JH+U));
Kc[eV0_,JH_,U_,t_]:=(2JH )/9 ( (t[[1]]-t[[3]])^2-3 t[[2]]^2 ) (eV0^2+3 JH^2-4 JH U+U^2)/((eV0+JH-U) (eV0+3 JH-U) (eV0-3 JH+U) (eV0-JH+U));
\[CapitalGamma]c[eV0_,JH_,U_,t_]:= (4 JH t[[2]] (t[[1]]-t[[3]]) )/9  (eV0^2+3 JH^2-4 JH U+U^2)/( (eV0+JH-U) (eV0+3 JH-U) (eV0-3 JH+U) (eV0-JH+U));
Jr[eV0_,JH_,U_,t_]:=Jc[eV0,JH,U,t]/Abs[Kc[0,JH,U,t]];
Kr[eV0_,JH_,U_,t_]:=Kc[eV0,JH,U,t]/Abs[Kc[0,JH,U,t]];
\[CapitalGamma]r[eV0_,JH_,U_,t_]:= \[CapitalGamma]c[eV0,JH,U,t]/Abs[Kc[0,JH,U,t]];
DMr[eV0_,JH_,U_,t_,Dmax_,s_:.8]:=Dmax eV0/ Abs[s(3 JH-U)]  Kc[eV0,JH,U,t]/Abs[Kc[s eV0,JH,U,t]];
JmatMicro[eV0_,JH_,U_,t_,Dmax_:.5,s_:.8]:=Jmat[Jr[eV0,JH,U,t]{1,1,1},Kr[eV0,JH,U,t]{1,1,1},\[CapitalGamma]r[eV0,JH,U,t]{1,1,1},{0,0,0},DMr[eV0,JH,U,t,Dmax,s]{1,1,1}];

addVortex[ Ko_,K0_, R_,L1_,L2_]:= Module[ {m=R[[1]],n=R[[2]],K=Ko,r=Table[0,6]},
K[[ m+n L1 +1,3    ]]=K0;  
K[[ Mod[m+1,L1]+Mod[n-1,L2] L1 +1  ,2    ]]=K0;
K[[ Mod[m+1,L1]+n L1 +1,1    ]]=K0;
K[[ Mod[m+1,L1]+Mod[n+1,L2] L1 +1  ,3   ]]=K0;
K[[ m+Mod[n+1,L2] L1+1   ,2    ]]=K0;
K[[ Mod[m-1,L1]+Mod[n+1,L2] L1 +1  ,1   ]]=K0;
K];
add4Vortices[ K0_,Kmod_, R1_,R2_, R3_,R4_,L_] :=Module[ {K}, K = addVortex[ K0,Kmod, R1,L,L]; K = addVortex[ K,Kmod, R2,L,L]; K = addVortex[ K,Kmod, R3,L,L]; addVortex[ K,Kmod, R4,L,L]];
Jarray4v[Jarray_,JmatMicro_,L_]:=Module[ {RS,RE,RW,RN,r2},   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor];
		RS={(r2-1),(r2-1)};
		RW={(r2-1),(L-r2-1)};
		RE={(L-r2-1),(r2-1)};
		RN={(L-r2-1),(L-r2-1)};
add4Vortices[ Jarray,JmatMicro, RS,RW, RE,RN,L] ];



(* ::Section::Bold:: *)
(*MF Loop -- microscopic parameters*)


ts = Table[ {5x,160,-12x,0,-60},{x,tV}];
hs =Table[  h[[1]]  hAngle[h[[2]],h[[3]]] , {h,hV}];  
eV0=0;U=2600;JH=300; dmax=.5; s0=1;

parametersMat=Table[Flatten[ Table[{
JmatMicro[     0,JH,U,ts[[t]],dmax,s0],hs[[h]],
JmatMicro[eV0,JH,U,ts[[t]],dmax,s0],hV[[h]],tV[[t]],eVs[[ev]]
} , {t,1,Length@tV},  {h,1,Length@hV}],1] ,  {ev,1,Length@eVs} ];



\[CapitalGamma]s=Table[x,{x,-.5,.5,.5}];
hV=Table[{x,0,0},{x,0.0001,1,.4}];
hs =Table[  h[[1]]  hAngle[h[[2]],h[[3]]] , {h,hV}];  
parametersMat={
Flatten[Table[{Jmat[0{1,1,1},-1{1,1,1},\[CapitalGamma]s[[t]]{1,1,1},0{1,1,1},0{1,1,1}], hs[[h]], {}, hV[[h]], 0, 0},{t,1,Length@\[CapitalGamma]s},{h,1,Length@hV}],1]
};


Print[" "];
Print["    NbName=",NbName,"; "];
Print["    Lambda s=", \[Eta]s,"; "];
Print["    Ls=",Ls,"; "];
Print["    tV=",tV,"; "];
Print["    ts=",ts,"; "];
Print["    hV=",hV,";"]
Print["    Couplings="];
Do[Print["    ", round[ fromJmat@parametersMat[[1,i,1]] ] ] , {i,1,Length@parametersMat[[1]] }];
(*Print["    Gammas=", \[CapitalGamma]s,"; "];*)
Print["    Steps=",steps,"; "];
Print["    acuracy=",acuracy,"; "];
(*
Print["    eVs=",eVs,"; "];
Print[" "];
Print["    Parameters="];
Do[
Print[ parameters[[i,j]] ],
{i,1,Length@parameters},
{j,1,Length@parameters[[i]] }
];*)
Print[" "];
tvf=AbsoluteTime[];
Print["Definition timing= ",round[tvf-t0] ]; 
t0=tvf;



(* ::Subsubsection::Bold::Closed:: *)
(*vortex free - Lagrange multiplier*)





(*Print[" "];Print[" "];Print["    Starting Lagrange multiplier free loop"];Print[" "];
t0=AbsoluteTime[]; *)


(*Do[
\[CapitalGamma]0=fromJmat[parametersMat[[1,p]][[1]]][[3]];
\[Alpha]0=-0.199+1.13 \[CapitalGamma]0;
\[Eta]s=Table[.1 \[Eta]+\[Alpha]0,{\[Eta],-1,1,.5}];
\[CapitalDelta]V\[Lambda]={};\[CapitalDelta]en\[Lambda]={};\[CapitalDelta]enSum\[Lambda]={};
Do[
Module[{ L,Jmat,Nc,h,\[CapitalLambda],T,H,\[Xi],EnG0,En,EnList={{},{},{}},u,u2,U,V,j,\[CapitalDelta]1=1,\[CapitalDelta]2=1,ES,gap,\[CapitalDelta]t,\[CapitalDelta]tHours,\[CapitalDelta]tMin,\[CapitalDelta]tSec,kTable,\[CapitalDelta]V,\[CapitalDelta]Vseq={},\[CapitalDelta]seq={},EMF,Esum,cMF,\[Eta],hp ,Tk=1/Sqrt[2] Tkmom,Umatvec,l=1}, L=Ls[[l]]; hp=Mod[p,Length@hV,1];
\[Eta]=\[Eta]s[[\[Eta]0]];
Print["Eta=",\[Eta]];
{Jmat,h}=parametersMat[[1,p]][[1;;2]];Nc=L^2;(*If[ p==1, UG=Uguess; VG=Vguess; ];*)  
UG=Uguess; VG=Vguess;  U=UG; V=VG;   kTable=toMomentumTable[L];

For[j=1,( ( j<(steps))\[And](Chop[ \[CapitalDelta]1,10^(-acuracy) ]!= 0) ), j++,  

	Umatvec=UmatVec[Jmat,h,U,V, kTable,Tk,\[Eta]] ;
u=Chop@Sum [  Module[{k,\[DoubleStruckCapitalU],\[DoubleStruckCapitalU]less,\[DoubleStruckCapitalU]gtr},k=kTable[[l]];\[DoubleStruckCapitalU]=Tk . Umatvec[[l]];\[DoubleStruckCapitalU]less=\[DoubleStruckCapitalU][[;;,5;;8]]; \[DoubleStruckCapitalU]gtr=\[DoubleStruckCapitalU][[;;,1;;4]];
(2I/Nc)Conjugate@Chop@{(*Conjugate[\[DoubleStruckCapitalU]gtr].Transpose[\[DoubleStruckCapitalU]gtr] Exp[I k.nx]+*)\[DoubleStruckCapitalU]less . \[DoubleStruckCapitalU]less\[ConjugateTranspose] Exp[-I k . nx], 
(*Conjugate[\[DoubleStruckCapitalU]gtr].Transpose[\[DoubleStruckCapitalU]gtr] Exp[I k.ny]+*)\[DoubleStruckCapitalU]less . \[DoubleStruckCapitalU]less\[ConjugateTranspose] Exp[-I k . ny], 
(*Conjugate[\[DoubleStruckCapitalU]gtr].Transpose[\[DoubleStruckCapitalU]gtr]+ *)\[DoubleStruckCapitalU]less . \[DoubleStruckCapitalU]less\[ConjugateTranspose]} 
  ],{l,1,Nc} ]; 

U[[1]]=u[[1]][[1;;4,5;;8]] ;U[[2]]=u[[2]][[1;;4,5;;8]];U[[3]]=u[[3]][[1;;4,5;;8]];V[[1]]=u[[3]][[1;;4,1;;4]];V[[2]]=u[[3]][[5;;8,5;;8]];
\[CapitalDelta]V=1/12 Sum[Sqrt@ Total[Power[#,2]&/@traceG[V[[\[Sigma]]]]  ],{\[Sigma],1,2}]; \[CapitalDelta]Vseq={\[CapitalDelta]Vseq,{j,\[CapitalDelta]V}} ;   
If[j>=2,\[CapitalDelta]2=\[CapitalDelta]1; \[CapitalDelta]1=Max[ Abs@(u-u2) ];\[CapitalDelta]seq={\[CapitalDelta]seq,{j,\[CapitalDelta]1}}  ;   ]; u2=u;                    
Print["    j=",j,"; \[CapitalDelta]V=",\[CapitalDelta]V,"; \[CapitalDelta]=",\[CapitalDelta]1];
                             ];
t1=AbsoluteTime[]; \[CapitalDelta]t=UnitConvert[ Quantity[t1 -t0, "Seconds" ], "Minutes" ]; t0=t1;
\[CapitalDelta]seq=Partition[Flatten[\[CapitalDelta]seq],2]; \[CapitalDelta]Vseq=Partition[Flatten[\[CapitalDelta]Vseq],2];
EMF   = enMFmom[Jmat,U,V,h,\[Eta]];
Esum  = enSUMmom[Jmat,U,V,h,L,\[Eta]];
cMF   = cMFmom[Jmat,U,V];
Print[];
Print[fromJmat[parametersMat[[1,p]][[1]]]];
Print["Eta=",\[Eta],(*"p=",p,"/",Length@parametersMat[[1]], *)"; j= ", j,"; Delta=",round@\[CapitalDelta]1,"; EMF=",EMF,"; Esum=",Esum,"; EMF-Esum=",round[EMF-Esum]   ];
Print[];
	(*dataToFile800[parametersMat[[1,p]],L,acuracy,{j,L,U,V,{0,0},{{EMF},{Esum},{cMF},\[CapitalDelta]seq,\[CapitalDelta]Vseq}},"free",NbName]; 

{jG,LG,\[Chi]G,\[Omega]G,\[Xi]G,EnG}= loadData[toPath800[parametersMat[[1,p]],L,acuracy,"free",NbName ]  ];*)
AppendTo[\[CapitalDelta]V\[Lambda],{\[Eta],\[CapitalDelta]V}];AppendTo[\[CapitalDelta]en\[Lambda],{\[Eta],EMF}];AppendTo[\[CapitalDelta]enSum\[Lambda],{\[Eta],Esum}];

 ]; , {\[Eta]0,1,Length@\[Eta]s}] ;
 
  dataToFile800[parametersMat[[1,p]],Ls[[1]],acuracy,
{fromJmat[parametersMat[[1,p]][[1]]],Ls[[1]],{0},{0},{0,0},{{0},{0},{0},{0},{0},{\[CapitalDelta]V\[Lambda],\[CapitalDelta]en\[Lambda],\[CapitalDelta]enSum\[Lambda],\[Eta]s}  }},
"free",NbName];  
 
 ,{p,1,Length@parametersMat[[1]] }] 
*)


(*dataToFile800[parametersMat[[1,1]],Ls[[1]],acuracy,
{0,Ls[[1]],{0},{0},{0,0},{{0},{0},{0},{0},{0},{\[CapitalDelta]V\[Lambda],\[CapitalDelta]en\[Lambda],\[CapitalDelta]enSum\[Lambda],\[Eta]s}  }},
"free",NbName]; *)



(*Module[{\[CapitalDelta]t}, t0v=AbsoluteTime[];
\[CapitalDelta]t= UnitConvert[ Quantity[N[t0v-tvf], "Seconds" ], "Minutes" ];   
Print["Free loop timing= ", IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]    ];t0=t0v;
Print[" "] Print[" "];];
Print["    Starting vortex free + electric field loop: "];Print[" "]*)


(* ::Subsubsection::Bold:: *)
(*vortex free*)


Print[" "];Print[" "];Print["    Starting free loop"];Print[" "];
t0=AbsoluteTime[]; 

Do[ \[CapitalGamma]0=fromJmat[parametersMat[[1,p]][[1]]][[3]];  \[Alpha]0=-0.2+1.13 \[CapitalGamma]0;
Module[{ L,Jmat,Nc,h,\[CapitalLambda],T,H,\[Xi],EnG0,En,EnList={{},{},{}},u,u2,U,V,j,\[CapitalDelta]1=1,\[CapitalDelta]2=1,ES,gap,\[CapitalDelta]t,\[CapitalDelta]tHours,\[CapitalDelta]tMin,\[CapitalDelta]tSec,kTable,\[CapitalDelta]V,\[CapitalDelta]Vseq={},\[CapitalDelta]seq={},EMF,Esum,cMF,\[Eta],hp ,Tk=1/Sqrt[2] Tkmom,Umatvec,l=1}, 
L=Ls[[l]]; hp=Mod[p,Length@hV,1];
\[Eta]=\[Alpha]0; Print["Eta=",\[Eta]];
{Jmat,h}=parametersMat[[1,p]][[1;;2]]; Nc=L^2;(*If[ p==1, UG=Uguess; VG=Vguess; ];*)  
UG=Uguess;  VG=Vguess;  U=UG;  V=VG;  kTable=toMomentumTable[L];

For[j=1,( (j<(steps))\[And](Chop[\[CapitalDelta]1,10^(-acuracy)]!=0) ), j++, Umatvec=UmatVec[Jmat,h,U,V,kTable,Tk,\[Eta]];
u=Chop@Sum[Module[{k,\[DoubleStruckCapitalU],\[DoubleStruckCapitalU]less,\[DoubleStruckCapitalU]gtr},k=kTable[[l]];\[DoubleStruckCapitalU]=Tk . Umatvec[[l]];\[DoubleStruckCapitalU]less=\[DoubleStruckCapitalU][[;;,5;;8]]; (*\[DoubleStruckCapitalU]gtr=\[DoubleStruckCapitalU][[;;,1;;4]];*)(2I/Nc)Conjugate@Chop@{\[DoubleStruckCapitalU]less . \[DoubleStruckCapitalU]less\[ConjugateTranspose] Exp[-I k . nx], \[DoubleStruckCapitalU]less . \[DoubleStruckCapitalU]less\[ConjugateTranspose] Exp[-I k . ny],\[DoubleStruckCapitalU]less . \[DoubleStruckCapitalU]less\[ConjugateTranspose]} ],{l,1,Nc} ]; 
U[[1]]=u[[1]][[1;;4,5;;8]]; U[[2]]=u[[2]][[1;;4,5;;8]]; U[[3]]=u[[3]][[1;;4,5;;8]]; V[[1]]=u[[3]][[1;;4,1;;4]]; V[[2]]=u[[3]][[5;;8,5;;8]];
\[CapitalDelta]V=1/(8Sqrt[3]) Sum[Sqrt@Total[  Power[#,2]&/@traceG[V[[\[Sigma]]]]  ],{\[Sigma],1,2}]; \[CapitalDelta]Vseq={\[CapitalDelta]Vseq,{j,\[CapitalDelta]V}} ;   
If[j>=2,\[CapitalDelta]2=\[CapitalDelta]1; \[CapitalDelta]1=Max[ Abs@(u-u2) ];\[CapitalDelta]seq={\[CapitalDelta]seq,{j,\[CapitalDelta]1}};   ];  u2=u;                    
Print["    j=",j,"; \[CapitalDelta]V=",\[CapitalDelta]V,"; \[CapitalDelta]=",\[CapitalDelta]1];                          ];
t1=AbsoluteTime[]; \[CapitalDelta]t=UnitConvert[ Quantity[t1-t0, "Seconds" ], "Minutes" ]; t0=t1;
\[CapitalDelta]seq=Partition[Flatten[\[CapitalDelta]seq],2]; \[CapitalDelta]Vseq=Partition[Flatten[\[CapitalDelta]Vseq],2];
EMF   = enMFmom[Jmat,U,V,h,\[Eta]];
Esum  = enSUMmom[Jmat,U,V,h,L,\[Eta]];
cMF   = cMFmom[Jmat,U,V];
Print[];
Print[fromJmat[parametersMat[[1,p]][[1]]]];
Print["Eta=",\[Eta],"p=",p,"/",Length@parametersMat[[1]],"; j= ", j,"; Delta=",round@\[CapitalDelta]1,"; EMF=",EMF,"; Esum=",Esum,"; EMF-Esum=",round[EMF-Esum]   ];
Print[];
	
dataToFile800[parametersMat[[1,p]],L,acuracy,{j,L,U,V,{0,0},{{EMF},{Esum},{cMF},\[CapitalDelta]seq,\[CapitalDelta]Vseq}},"free",NbName]; 
(*{jG,LG,\[Chi]G,\[Omega]G,\[Xi]G,EnG}= loadData[toPath800[parametersMat[[1,p]],L,acuracy,"free",NbName ]  ];  AppendTo[\[CapitalDelta]V\[Lambda],{\[Eta],\[CapitalDelta]V}];AppendTo[\[CapitalDelta]en\[Lambda],{\[Eta],EMF}];AppendTo[\[CapitalDelta]enSum\[Lambda],{\[Eta],Esum}]; *)  ];   
,{p,1,Length@parametersMat[[1]]} ]; Print[" "]; 



Module[{\[CapitalDelta]t}, t0v=AbsoluteTime[];
\[CapitalDelta]t= UnitConvert[ Quantity[N[t0v-tvf], "Seconds" ], "Minutes" ];   
Print["Free loop timing= ", IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]    ];t0=t0v;
Print[" "];
Print[" "];   ];
Print["    Starting vortex free + electric field loop: "];Print[" "]


(* ::Subsubsection::Bold:: *)
(*vortex free + gradually increase  electric field*)


(*Module[{\[CapitalDelta]t}, t0v=AbsoluteTime[];
\[CapitalDelta]t= UnitConvert[ Quantity[N[t0v-tvf], "Seconds" ], "Minutes" ];   
Print["Free loop timing= ", IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]    ];t0=t0v;
Print[" "] Print[" "];];
Print["    Starting vortex free + electric field loop: "];Print[" "]*)


(* ::Subsubsection::Bold:: *)
(*four vortex + gradually increase  electric field*)


(*Module[{\[CapitalDelta]t},t4v=AbsoluteTime[];\[CapitalDelta]t= UnitConvert[ Quantity[N[t4v-t0v], "Seconds" ], "Hours" ];
\[CapitalDelta]tHours=IntegerPart[\[CapitalDelta]t];
\[CapitalDelta]tMin=IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ];
\[CapitalDelta]tSec=IntegerPart@UnitConvert[FractionalPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ], "Seconds" ];
Print[ "Free loop + electric field timing= ",ToString@\[CapitalDelta]tHours," : ",ToString@\[CapitalDelta]tMin," : ",ToString@\[CapitalDelta]tSec    ];
t0=t4v;Print[" "] Print[" "];
];
Print["    Starting four vortex + electric field loop: "];Print[" "]*)


(*Module[{\[CapitalDelta]t},t1=AbsoluteTime[];\[CapitalDelta]t= UnitConvert[ Quantity[N[t1-t4v], "Seconds" ], "Hours" ];
\[CapitalDelta]tHours=IntegerPart[\[CapitalDelta]t];
\[CapitalDelta]tMin=IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ];
\[CapitalDelta]tSec=IntegerPart@UnitConvert[FractionalPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ], "Seconds" ];
Print[ "4 vortices loop timing \[CapitalDelta]t = ",ToString@\[CapitalDelta]tHours," : ",ToString@\[CapitalDelta]tMin," : ",ToString@\[CapitalDelta]tSec    ] 
];*)


(* ::Subsubsection::Bold:: *)
(*vortex free + gradually increase  parameters *)


(*Module[{\[CapitalDelta]t}, t0v=AbsoluteTime[];
\[CapitalDelta]t= UnitConvert[ Quantity[N[t0v-tvf], "Seconds" ], "Minutes" ];   
Print["Free loop timing= ", IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]    ];t0=t0v;Print[" "] Print[" "];];
Print["    Starting vortex free + electric field loop: "];Print[" "]*)


(* ::Subsubsection::Bold:: *)
(*four vortex + gradually increase  parameters *)


(*Module[{\[CapitalDelta]t},t4v=AbsoluteTime[];
\[CapitalDelta]t= UnitConvert[ Quantity[N[t4v-t0v], "Seconds" ], "Minutes" ];   
Print["Free loop + electric field timing= ", IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]    ];t0=t4v;Print[" "] Print[" "];];
Print["    Starting four vortex -fixed eV- gradually changing parameters  "];Print[" "]*)



CloseKernels[];


(*Module[{\[CapitalDelta]t},t1=AbsoluteTime[];\[CapitalDelta]t= UnitConvert[ Quantity[N[t1-t4v], "Seconds" ], "Hours" ];
\[CapitalDelta]tHours=IntegerPart[\[CapitalDelta]t];
\[CapitalDelta]tMin=IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ];
\[CapitalDelta]tSec=IntegerPart@UnitConvert[FractionalPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ], "Seconds" ];
Print[ "4 vortices loop timing \[CapitalDelta]t = ",ToString@\[CapitalDelta]tHours," : ",ToString@\[CapitalDelta]tMin," : ",ToString@\[CapitalDelta]tSec    ] ];*)
