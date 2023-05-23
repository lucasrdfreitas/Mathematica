(* ::Package:: *)

t0=AbsoluteTime[];
Print["Before Starting Kernels"];


Needs["ClusterIntegration`"];
(*kernels = aunchKernels[SGE["micro4", 10]];*)
Quiet[kernels = LaunchKernels[]];
Print["Starting Kernels"];


(* ::Section::Bold:: *)
(*Definitions *)


(* ::Code::Bold:: *)
(**)


NbName="705"; \[Lambda]0=0.5; 

		Ls = Range[32,32,4]; 				tV={3};				
		hV={ (*{0.25,0,90},{0.25,15,90},{0.25,30,90},{0.25,45,90},{0.25,60,90},*){0.25,75,90},{0.25,90,90} 
		(*,{0.2612,45,45},{0.2612,45,90},{0.2612,90,0},{0.2612,90,45}*)   };

		steps=500;				acuracy=6;     eVs=Table[1700 x, {x,0,0,0.0499999}];  (* eV=\[Xi](U-3JH)=1500\[Xi] *)


(* ::Subsubsection::Bold::Closed:: *)
(*Constants*)


(* ::Code::Bold:: *)
(**)


one[i_,j_,N_]:= SparseArray[ {i,j} -> 1, {N,N}];  block[i_,j_] := one[i,j,8];

hc=1/Sqrt[3] {1,1,1};hb=1/Sqrt[2] {1,-1,0};ha=1/Sqrt[6] {1,1,-2};hx={1,0,0};
hAngle[\[Theta]_,\[Phi]_]:= ha Cos[\[Phi] \[Pi]/180] Sin[\[Theta] \[Pi]/180]+hb Sin[\[Phi] \[Pi]/180] Sin[\[Theta] \[Pi]/180]+hc Cos[\[Theta] \[Pi]/180];  (* Magnetic field directions a, b, c*)
nx={1/2,Sqrt[3]/2};  ny={-(1/2),Sqrt[3]/2};                                                            (* Honeycomb Bravais basis  *)

ReverseSort[list_]:=Reverse@Sort@list;
round[\[Kappa]_]:=N[Round[10000\[Kappa]]/10000];round\[CapitalDelta][\[Kappa]_]:=N[Round[10^13 \[Kappa]]/10^13];

\[Delta]x={-Sqrt[3]/2,-1/2}/Sqrt[3];      
\[Delta]y={Sqrt[3]/2,-1/2}/Sqrt[3];     
\[Delta]z={0,1}/Sqrt[3];   
(* NN vectors  *)
\[Chi]Gx = {{0.5249,0,0,0},{0,-1,0,0},{0,0,0,0},{0,0,0,0}};\[Chi]Gy = {{0.5249,0,0,0},{0,0,0,0},{0,0,-1,0},{0,0,0,0}};\[Chi]Gz ={{0.5249,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,-1}};

\[Omega]GA = {{I,0.000001,0.000001,0.000001},{-0.000001,I,0.000001,-0.000001},{-0.000001,-0.000001,I,0.000001},{-0.000001,0.000001,-0.000001,I}}; \[Omega]GB = \[Omega]GA;


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
createDir@FileNameJoin[{Directory[],"Files","pure",gauge, StringReplace["t=X1_eV=X2_JKG=X3_JKGmod=X4",
{"X1"->ToString[parameters[[9]]],"X2"->ToString[parameters[[10]]],"X3"->ToString[parameters[[1;;3]]],"X4"->ToString[parameters[[5;;7]]] }] }];
		path = toPathPure[parameters,L,acuracy,gauge];		
		(*Print["Pure path=",path];Print[];*)
		f = OpenAppend[path];
		 Write[ f, data];
		 Close[f];                ];
		 
		 
toPathPure[parameters0_,L_,acuracy_,gauge_]:= Module[{h ,hS,parameters=parameters0, r,\[Phi],\[Theta]},h=parameters[[4]] ;hS=parameters[[8]]; {r,\[Theta],\[Phi]}=hS;   
parameters[[1]]=0;parameters[[3]]=0;parameters[[5]]=0;parameters[[7]]=0;
FileNameJoin[{Directory[], "Files" ,"pure",gauge,StringReplace["t=X1_eV=X2_JKG=X3_JKGmod=X4",
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


(* ::Subsection::Bold::Closed:: *)
(*MF definitions*)


(* ::Subsubsection::Bold::Closed:: *)
(*Saving and Loading data*)


createDir[path_] :=
Module[ {l=Length@FileNames[path]},		
	If[ l==0,  (*CreateDirectory@File@FileNameJoin[{path}];*) 
	CreateDirectory@File@FileNameJoin[{path,"data" }]; CreateDirectory@File@FileNameJoin[{path,"graph" }]; , Null]  ];

loadData[pathData_]:=
Module[ {f,data},		        f = OpenRead[pathData];
If[f==$Failed, Print["Failed to OpenRead file at: ", pathData ];Return[$Failed]; Abort[] ];
		        data=ReadList[f];
		        Close[f];			data[[-1]]
];

loadDataTry[pathData_]:=
Module[ {f,data,path=FindFile[pathData]},
If[path==$Failed,(*Print["New entry at:",pathData];*)Return[$Failed]];
		        f = OpenRead[pathData];
(*If[f==$Failed, Print["Failed to OpenRead file at: ", pathData ]; Abort[] ];*)
		        data=ReadList[f];
		        Close[f];			data[[-1]]
];

toPath [parameters0_,L_,acuracy_,gauge_,NbName_]:= Module[{h ,hS,parameters=parameters0, r,\[Phi],\[Theta]},h=parameters[[4]] ;hS=parameters[[8]]; {r,\[Theta],\[Phi]}=hS;   
If[gauge=="free",parameters[[5;;7]]=parameters0[[1;;3]];parameters[[10]]=0;];
FileNameJoin[{    Directory[]   , "Files" ,NbName,gauge,StringReplace["t=X1_eV=X2_JKG=X3_JKGmod=X4",{"X1"->  ToString[parameters[[9]]],"X2"->  ToString[parameters[[10]]],"X3"->  ToString[parameters[[1;;3]]  ],"X4"->  ToString[parameters[[5;;7]]]   }] , "data"  , 
StringReplace["h=(M,N,T)_L=Y_A=Z.txt",{"Y"-> ToString[L], "Z"-> ToString[acuracy],"M"->  ToString[r,InputForm] ,"N"->ToString@\[Phi],"T"->ToString@\[Theta]}   ] 
   }]];
   
dataToFile[ parameters0_,L_,acuracy_,data_,gauge_,NbName_] :=
Module[ {path,f,parameters=parameters0}, If[gauge=="free",parameters[[5;;7]]=parameters0[[1;;3]];parameters[[10]]=0;];
		createDir@FileNameJoin[{ Directory[] ,"Files",  NbName,gauge, StringReplace["t=X1_eV=X2_JKG=X3_JKGmod=X4",
{"X1"->  ToString[parameters[[9]]],"X2"->  ToString[parameters[[10]]],"X3"->  ToString[parameters[[1;;3]]  ],"X4"->  ToString[parameters[[5;;7]]]   }]      }] ;
		path = toPath[ parameters,L,acuracy,gauge,NbName];		
		f = OpenAppend[path];
		 Write[ f, data];
		 Close[f];                ];


(* ::Subsubsection::Bold::Closed:: *)
(*Hamiltonian  matrices*)


(* ::Code::Bold:: *)
(**)


auxHx[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_,r_]:={
{0,0,0,0,-K[[r,1]] \[Chi][[1,r,2,2]],K[[r,1]] \[Chi][[1,r, 2,1]],0,0},
{0,0,0,0, K[[r,1]] \[Chi][[1,r, 1,2]],-K[[r,1]] \[Chi][[1,r, 1,1]],0,0},
{0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,J[[r,1]] (-\[Chi][[1, r,2,2]]-\[Chi][[1, r,3,3]]-\[Chi][[1,r, 4,4]]),J[[r,1]] \[Chi][[1,r, 2,1]],J[[r,1]] \[Chi][[1,r, 3,1]],J[[r,1]] \[Chi][[1,r, 4,1]]},
{0,0,0,0,J[[r,1]] \[Chi][[1,r, 1,2]],-J[[r,1]] \[Chi][[1,r, 1,1]],0,0},
{0,0,0,0,J[[r,1]] \[Chi][[1,r, 1,3]],0,-J[[r,1]] \[Chi][[1,r, 1,1]],0},
{0,0,0,0,J[[r,1]] \[Chi][[1,r, 1,4]],0,0,-J[[r,1]] \[Chi][[1,r, 1,1]]},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,\[CapitalGamma][[r,1]] (-\[Chi][[1,r,3,4]]-\[Chi][[1,r, 4,3]]),0,\[CapitalGamma][[r,1]] \[Chi][[1,r, 4,1]],\[CapitalGamma][[r,1]]  \[Chi][[1,r, 3,1]]},
{0,0,0,0,0,0,0,0},
{0,0,0,0,\[CapitalGamma][[r,1]] \[Chi][[1,r, 1,4]],0,0,\[CapitalGamma][[r,1]] \[Chi][[1,r, 1,1]]},
{0,0,0,0,\[CapitalGamma][[r,1]] \[Chi][[1,r, 1,3]],0,\[CapitalGamma][[r,1]] \[Chi][[1,r, 1,1]],0},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}; 

auxHy[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_,r_]:={
{0,0,0,0,-K[[r,2]] \[Chi][[2,r, 3,3]],0,K [[r,2]] \[Chi][[2,r, 3,1]],0},
{0,0,0,0,0,0,0,0},
{0,0,0,0,K[[r,2]] \[Chi][[2,r, 1,3]],0,-K [[r,2]] \[Chi][[2,r, 1,1]],0},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,J[[r,2]](-\[Chi][[2,r, 2,2]]-\[Chi][[2,r, 3,3]]-\[Chi][[2,r, 4,4]]),J[[r,2]]  \[Chi][[2,r, 2,1]],J[[r,2]]  \[Chi][[2,r, 3,1]],J[[r,2]]  \[Chi][[2,r, 4,1]]},
{0,0,0,0,J[[r,2]]\[Chi][[2,r, 1,2]],-J[[r,2]] \[Chi][[2,r, 1,1]],0,0},
{0,0,0,0,J[[r,2]] \[Chi][[2,r, 1,3]],0,-J[[r,2]] \[Chi][[2,r, 1,1]],0},
{0,0,0,0,J[[r,2]] \[Chi][[2,r, 1,4]],0,0,-J[[r,2]] \[Chi][[2,r, 1,1]]},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,\[CapitalGamma][[r,2]] (-\[Chi][[2,r, 2,4]]-\[Chi][[2, r,4,2]]),\[CapitalGamma][[r,2]] \[Chi][[2,r, 4,1]],0,\[CapitalGamma][[r,2]] \[Chi][[2,r, 2,1]]},
{0,0,0,0,\[CapitalGamma][[r,2]] \[Chi][[2,r, 1,4]],0,0,\[CapitalGamma][[r,2]] \[Chi][[2,r, 1,1]]},
{0,0,0,0,0,0,0,0},
{0,0,0,0,\[CapitalGamma][[r,2]] \[Chi][[2,r, 1,2]],\[CapitalGamma][[r,2]] \[Chi][[2,r, 1,1]],0,0},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}};

HeffList[J_,K_,\[CapitalGamma]_,h_,\[Omega]_]:=Module[{L,Nc},Nc=Length@K; L= Round[Sqrt[Nc] ];      Table[
Module[ {RA,RB,m,n,mAx,nAy,mBx,nBy},  n=\[LeftFloor](r-1)/L\[RightFloor];m=r-1-n L;  
mAx=Mod[m+1,L];nAy=Mod[n+1,L]; mBx=Mod[m-1,L]; nBy=Mod[n-1,L];  RA={mAx+n L+1,m+nAy L+1,m+n L+1}; RB={mBx+n L+1,m+nBy L+1,m+n L+1}; 

{ {   h[[1]]-K [[r,1]] \[Omega][[2,r,2,1]]-(J[[r,1]]+J[[r,2]]+J[[r,3]]) \[Omega][[2,r,2,1]] -2\[CapitalGamma][[r,2]] \[Omega][[2,RA[[2]],4,1]] -2\[CapitalGamma][[r,3]] \[Omega][[2,RA[[3]],3,1]] 
,h[[2]]-K [[r,2]] \[Omega][[2,r,3,1]]-(J[[r,1]]+J[[r,2]]+J[[r,3]]) \[Omega][[2,r,3,1]] -2\[CapitalGamma][[r,3]] \[Omega][[2,RA[[3]],2,1]] -2\[CapitalGamma][[r,1]] \[Omega][[2,RA[[1]],4,1]]  
,h[[3]]-K [[r,3]] \[Omega][[2,r,4,1]]-(J[[r,1]]+J[[r,2]]+J[[r,3]]) \[Omega][[2,r,4,1]] -2\[CapitalGamma][[r,1]] \[Omega][[2,RA[[1]],3,1]] -2\[CapitalGamma][[r,2]] \[Omega][[2,RA[[2]],2,1]] 
} ,
{    h[[1]]-K [[RB[[1]],1]] \[Omega][[1,RB[[1]],2,1]]-Sum[J[[RB[[\[Beta]]],\[Beta]]]\[Omega][[1,RB[[\[Beta]]],2,1]],{\[Beta],1,3}] -2\[CapitalGamma][[RB[[2]],2]] \[Omega][[1,RB[[2]],4,1]] -2\[CapitalGamma][[RB[[3]],3]] \[Omega][[1,RB[[3]],3,1]] 
,h[[2]]-K [[RB[[2]],2]] \[Omega][[1,RB[[2]],3,1]]-Sum[J[[RB[[\[Beta]]],\[Beta]]]\[Omega][[1,RB[[\[Beta]]],3,1]],{\[Beta],1,3}] -2\[CapitalGamma][[RB[[3]],3]] \[Omega][[1,RB[[3]],2,1]] -2\[CapitalGamma][[RB[[1]],1]] \[Omega][[1,RB[[3]],4,1]]  
,h[[3]]-K [[RB[[3]],3]] \[Omega][[1,RB[[3]],4,1]]-Sum[J[[RB[[\[Beta]]],\[Beta]]]\[Omega][[1,RB[[\[Beta]]],4,1]],{\[Beta],1,3}] -2\[CapitalGamma][[RB[[1]],1]] \[Omega][[1,RB[[1]],3,1]] -2\[CapitalGamma][[RB[[2]],2]] \[Omega][[1,RB[[2]],2,1]] 
}}
],    {r,1,Nc}]    ];

(*\[Lambda]2List[Heff_,\[Omega]_]:=Module[{L,Nc},Nc=Length@Heff; L= Round[Sqrt[Nc] ];
Table[Module[ {RA,RB,m,n,mAx,nAy,mBx,nBy},  n=\[LeftFloor](r-1)/L\[RightFloor];m=r-1-n L;  
mAx=Mod[m+1,L];nAy=Mod[n+1,L]; mBx=Mod[m-1,L];nBy=Mod[n-1,L];  RA={mAx+n L+1,m+nAy L+1,m+n L+1}; RB={mBx+n L+1,m+nBy L+1,m+n L+1}; 
{{   Heff[[r,1,1]] \[Omega][[2,r,2,1]] /(\[Omega][[2,r,2,1]] -\[Omega][[2,r,3,4]] ),   Heff[[r,1,2]] \[Omega][[2,r,3,1]] /(\[Omega][[2,r,3,1]] -\[Omega][[2,r,4,2]] ),   Heff[[r,1,3]] \[Omega][[2,r,4,1]] /(\[Omega][[2,r,4,1]] -\[Omega][[2,r,2,3]] )} ,
{   Heff[[r,2,1]] \[Omega][[1,r,2,1]] /(\[Omega][[1,r,2,1]] -\[Omega][[1,r,3,4]] ),   Heff[[r,2,2]] \[Omega][[1,r,3,1]] /(\[Omega][[1,r,3,1]] -\[Omega][[1,r,4,2]] ),   Heff[[r,2,3]] \[Omega][[1,r,4,1]] /(\[Omega][[1,r,4,1]] -\[Omega][[1,r,2,3]] )}}
],    {r,1,Nc}] ];
\[Lambda]1List[Heff_,\[Omega]_]:=Module[{L,Nc},Nc=Length@Heff; L= Round[Sqrt[Nc] ];
Table[
{{   Heff[[r,1,1]] (-\[Omega][[2,r,3,4]] )/(\[Omega][[2,r,2,1]] -\[Omega][[2,r,3,4]] ),   Heff[[r,1,2]] (-\[Omega][[2,r,4,2]])/(\[Omega][[2,r,3,1]] -\[Omega][[2,r,4,2]] ),   Heff[[r,1,3]] (-\[Omega][[2,r,2,3]] )/(\[Omega][[2,r,4,1]] -\[Omega][[2,r,2,3]] )} ,
{   Heff[[r,2,1]] (-\[Omega][[1,r,3,4]])/(\[Omega][[1,r,2,1]] -\[Omega][[1,r,3,4]] ),   Heff[[r,2,2]] (-\[Omega][[1,r,4,2]])/(\[Omega][[1,r,3,1]] -\[Omega][[1,r,4,2]] ),   Heff[[r,2,3]] (-\[Omega][[1,r,2,3]] )/(\[Omega][[1,r,4,1]] -\[Omega][[1,r,2,3]] )}}
,    {r,1,Nc}] ];*)
\[Lambda]2List[Heff_,\[Omega]_]:=Module[{L,Nc},Nc=Length@Heff; L= Round[Sqrt[Nc] ];
Table[Module[ {RA,RB,m,n,mAx,nAy,mBx,nBy},  n=\[LeftFloor](r-1)/L\[RightFloor];m=r-1-n L;  
mAx=Mod[m+1,L];nAy=Mod[n+1,L]; mBx=Mod[m-1,L];nBy=Mod[n-1,L];  RA={mAx+n L+1,m+nAy L+1,m+n L+1}; RB={mBx+n L+1,m+nBy L+1,m+n L+1}; 
Re@{{   Heff[[r,1,1]] \[Omega][[2,r,2,1]] /(\[Omega][[2,r,2,1]] -\[Omega][[2,r,3,4]]+I (10^-9.9) ),   Heff[[r,1,2]] \[Omega][[2,r,3,1]] /(\[Omega][[2,r,3,1]] -\[Omega][[2,r,4,2]]+I (10^-9.9) ),   Heff[[r,1,3]] \[Omega][[2,r,4,1]] /(\[Omega][[2,r,4,1]] -\[Omega][[2,r,2,3]] +I 10^-9.9)} ,
{   Heff[[r,2,1]] \[Omega][[1,r,2,1]] /(\[Omega][[1,r,2,1]] -\[Omega][[1,r,3,4]]+I (10^-9.9) ),   Heff[[r,2,2]] \[Omega][[1,r,3,1]] /(\[Omega][[1,r,3,1]] -\[Omega][[1,r,4,2]] +I 10^-9.9),   Heff[[r,2,3]] \[Omega][[1,r,4,1]] /(\[Omega][[1,r,4,1]] -\[Omega][[1,r,2,3]]+I 10^-9.9)}}
],    {r,1,Nc}] ];
\[Lambda]1List[Heff_,\[Omega]_]:=Module[{L,Nc},Nc=Length@Heff; L= Round[Sqrt[Nc] ];
Table[
Re@{{   Heff[[r,1,1]] (-\[Omega][[2,r,3,4]] )/(\[Omega][[2,r,2,1]] -\[Omega][[2,r,3,4]]+I (10^-9.9) ),   Heff[[r,1,2]] (-\[Omega][[2,r,4,2]])/(\[Omega][[2,r,3,1]] -\[Omega][[2,r,4,2]] +I 10^-9.9),   Heff[[r,1,3]] (-\[Omega][[2,r,2,3]] )/(\[Omega][[2,r,4,1]] -\[Omega][[2,r,2,3]]+I (10^-9.9) )} ,
{   Heff[[r,2,1]] (-\[Omega][[1,r,3,4]])/(\[Omega][[1,r,2,1]] -\[Omega][[1,r,3,4]]+I (10^-9.9) ),   Heff[[r,2,2]] (-\[Omega][[1,r,4,2]])/(\[Omega][[1,r,3,1]] -\[Omega][[1,r,4,2]] +I 10^-9.9),   Heff[[r,2,3]] (-\[Omega][[1,r,2,3]] )/(\[Omega][[1,r,4,1]] -\[Omega][[1,r,2,3]]+I (10^-9.9) )}}
,    {r,1,Nc}] ];

auxHz[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_,r_,\[Lambda]1_,\[Lambda]2_,heff_]:={
{0,0,0,0,         -K[[r,3]] \[Chi][[3,r, 4,4]],0,0, K[[r,3]] \[Chi][[3,r, 4,1]]    },
{0,0,0,0,         0,0,0,0   },
{0,0,0,0,         0,0,0,0   },
{0,0,0,0,         K[[r,3]] \[Chi][[3,r, 1,4]],0,0,-K[[r,3]] \[Chi][[3,r, 1,1]]  },
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,         J[[r,3]](-\[Chi][[3, r,2,2]]-\[Chi][[3,r, 3,3]]-\[Chi][[3,r, 4,4]]),J[[r,3]] \[Chi][[3, r,2,1]],J[[r,3]] \[Chi][[3,r, 3,1]],J[[r,3]] \[Chi][[3,r, 4,1]]    },
{0,0,0,0,         J[[r,3]] \[Chi][[3,r, 1,2]],-J[[r,3]] \[Chi][[3,r, 1,1]],0,0   },
{0,0,0,0,         J[[r,3]] \[Chi][[3,r, 1,3]],0,-J[[r,3]] \[Chi][[3,r, 1,1]], 0  },
{0,0,0,0,         J[[r,3]] \[Chi][[3,r, 1,4]],0,0,-J[[r,3]] \[Chi][[3,r, 1,1]]  },
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,         \[CapitalGamma][[r,3]](-\[Chi][[3,r, 2,3]]-\[Chi][[3,r, 3,2]]),\[CapitalGamma][[r,3]] \[Chi][[3,r, 3,1]],\[CapitalGamma][[r,3]] \[Chi][[3,r, 2,1]],0   },
{0,0,0,0,         \[CapitalGamma][[r,3]] \[Chi][[3,r, 1,3]],0,\[CapitalGamma][[r,3]] \[Chi][[3,r, 1,1]],0   },
{0,0,0,0,         \[CapitalGamma][[r,3]] \[Chi][[3,r, 1,2]],\[CapitalGamma][[r,3]] \[Chi][[3,r, 1,1]],0, 0  },
{0,0,0,0,         0,0,0,0  },
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,             0,0,0,0 },
{-heff[[r,2,2]]+\[Lambda]1[[r,2,1]] ,0,\[Lambda]2[[r,2,3]],0,         0,0,0,0   },
{-heff[[r,2,2]]+\[Lambda]1[[r,2,2]] ,0,0,\[Lambda]2[[r,2,1]],         0,0,0,0   },
{-heff[[r,2,3]]+\[Lambda]1[[r,2,3]] ,\[Lambda]2[[r,2,2]],0,0,         0,0,0,0   },
{0,0,0,0,0,0,0,0},
{0,0,0,0,-heff[[r,1,1]]+\[Lambda]1[[r,1,1]]  ,0,\[Lambda]2[[r,1,3]],0},
{0,0,0,0,-heff[[r,1,2]]+\[Lambda]1[[r,1,2]] ,0,0,\[Lambda]2[[r,1,1]]},
{0,0,0,0,-heff[[r,1,3]]+\[Lambda]1[[r,1,3]] ,\[Lambda]2[[r,1,2]],0,0} }; 
					
HMF[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_,L1_,L2_,\[Lambda]1_,\[Lambda]2_,Heff_] := Module[{Nc=L1 L2,H}, H= I Sum[Module[{rx,ry,rz,mx,ny}, mx=Mod[m+1,L1,1];ny=Mod[n+1,L2,1];
 rx=mx+(n-1)L1;ry=m+(ny-1)L1;rz=m+(n-1)L1;
 KroneckerProduct[ one[rz,rx,Nc],auxHx[J,K,\[CapitalGamma],h,\[Chi],\[Omega],rz]]+  KroneckerProduct[ one[rz,ry,Nc],auxHy[J,K,\[CapitalGamma],h,\[Chi],\[Omega],rz]]+ KroneckerProduct[ one[rz,rz,Nc],auxHz[J,K,\[CapitalGamma],h,\[Chi],\[Omega],rz,\[Lambda]1,\[Lambda]2,Heff] ]
] ,{m,1,L1},{n,1,L2}];  H+H\[ConjugateTranspose]];



(* ::Text::Bold:: *)
(*Energy*)


(* ::Code::Bold:: *)
(**)


EnMF[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_,L1_,L2_] := Module[{Nc=L1 L2,E},
E= Sum[Module[{r,mx,ny,\[CurlyEpsilon]=Abs@LeviCivitaTensor[3]}, 	mx=Mod[m+1,L1];ny=Mod[n+1,L2];   r={mx+n L1+1,m+ny L1+1,m+n L1+1};
-h[[ \[Gamma]]] (   \[Omega][[1,r[[3]], \[Gamma]+1,1]]+\[Omega][[2,r[[ 3]] , \[Gamma]+1,1]]    )         

+J[[r[[3]]  ,\[Gamma]]]  Sum[    \[Omega][[1,r[[3]],\[Beta]+1,1]]\[Omega][[2,r[[ \[Gamma]]] ,\[Beta]+1,1]]  - \[Chi][[ \[Gamma],r[[3]],\[Beta]+1,\[Beta]+1]]   \[Chi][[ \[Gamma],r[[3]],1,1]]  +\[Chi][[ \[Gamma],r[[3]],\[Beta]+1,1]]   \[Chi][[ \[Gamma],r[[3]],1,1+\[Beta]]]      ,{\[Beta],1,3}] 

+K[[r[[3]]  , \[Gamma]]] (  \[Omega][[1,r[[3]], \[Gamma]+1,1]]\[Omega][[2,r[[ \[Gamma]]] , \[Gamma]+1,1]]  - \[Chi][[ \[Gamma],r[[3]], \[Gamma]+1, \[Gamma]+1]]   \[Chi][[ \[Gamma],r[[3]],1,1]]  +\[Chi][[ \[Gamma],r[[3]], \[Gamma]+1,1]]   \[Chi][[ \[Gamma],r[[3]],1,1+ \[Gamma]]]    )

+Sum[    \[CapitalGamma][[r[[3]]  ,\[Gamma]]]  \[CurlyEpsilon][[\[Alpha],\[Beta],\[Gamma]]] (  \[Omega][[1,r[[3]], \[Alpha]+1,1]]\[Omega][[2,r[[ \[Gamma]]] , \[Beta]+1,1]]  - \[Chi][[ \[Gamma],r[[3]], \[Alpha]+1, \[Beta]+1]]   \[Chi][[ \[Gamma],r[[3]],1,1]]  +\[Chi][[ \[Gamma],r[[3]], \[Alpha]+1,1]]   \[Chi][[ \[Gamma],r[[3]],1,1+ \[Beta]]]    )
,     {\[Alpha],1,3},  {\[Beta],1,3}]    ]

,     { \[Gamma],1,3},{m,0,L1-1},{n,0,L2-1} ]; Chop@(E/ Nc )   ];



EnMF0[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_,L1_,L2_] := Module[{Nc=L1 L2,E},
E= Sum[Module[{r,mx,ny,\[CurlyEpsilon]=Abs@LeviCivitaTensor[3]}, 	mx=Mod[m+1,L1];ny=Mod[n+1,L2];   r={mx+n L1+1,m+ny L1+1,m+n L1+1};
-h[[ \[Gamma]]] (   \[Omega][[1,r[[3]], \[Gamma]+1,1]]+\[Omega][[2,r[[ 3]] , \[Gamma]+1,1]]    )         

+J[[r[[3]] ,\[Gamma]]]  Sum[    \[Omega][[1,r[[3]],\[Beta]+1,1]]\[Omega][[2,r[[ \[Gamma]]] ,\[Beta]+1,1]]  +Abs[- \[Chi][[ \[Gamma],r[[3]],\[Beta]+1,\[Beta]+1]]   \[Chi][[ \[Gamma],r[[3]],1,1]] ]  +\[Chi][[ \[Gamma],r[[3]],\[Beta]+1,1]]   \[Chi][[ \[Gamma],r[[3]],1,1+\[Beta]]]      ,{\[Beta],1,3}] 

+K[[r[[3]] ,\[Gamma]]] (  \[Omega][[1,r[[3]], \[Gamma]+1,1]]\[Omega][[2,r[[ \[Gamma]]] , \[Gamma]+1,1]]  + Abs[- \[Chi][[ \[Gamma],r[[3]], \[Gamma]+1, \[Gamma]+1]]   \[Chi][[ \[Gamma],r[[3]],1,1]] ] +\[Chi][[ \[Gamma],r[[3]], \[Gamma]+1,1]]   \[Chi][[ \[Gamma],r[[3]],1,1+ \[Gamma]]]    )

+Sum[    2\[CapitalGamma][[r[[3]] ,\[Gamma]]]  \[CurlyEpsilon][[\[Alpha],\[Beta],\[Gamma]]] (  \[Omega][[1,r[[3]], \[Alpha]+1,1]]\[Omega][[2,r[[ \[Gamma]]] , \[Beta]+1,1]]  +Abs[- \[Chi][[ \[Gamma],r[[3]], \[Alpha]+1, \[Beta]+1]]   \[Chi][[ \[Gamma],r[[3]],1,1]] ] +\[Chi][[ \[Gamma],r[[3]], \[Alpha]+1,1]]   \[Chi][[ \[Gamma],r[[3]],1,1+ \[Beta]]]    )
,     {\[Alpha],1,3},  {\[Beta],1,3}]    ]

,     { \[Gamma],1,3},{m,0,L1-1},{n,0,L2-1} ]; Chop@(E/ Nc )   ];

EnLagMF[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_,L1_,L2_,\[Lambda]_:\[Lambda]0] := Module[{Nc=L1 L2,E},
E= 1/2 \[Lambda] Sum[Module[{r,mx,ny}, 	mx=Mod[m+1,L1];ny=Mod[n+1,L2];   r={mx+n L1+1,m+ny L1+1,m+n L1+1};
      h[[1]] (   \[Omega][[\[Sigma],r[[3]],1+1,1]]+\[Omega][[\[Sigma],r[[3]] , 2+1,3+1]]    )  
+h[[2]] (   \[Omega][[\[Sigma],r[[3]],2+1,1]]+\[Omega][[\[Sigma],r[[3]] , 3+1,1+1]]    )  
+h[[3]] (   \[Omega][[\[Sigma],r[[3]],3+1,1]]+\[Omega][[\[Sigma],r[[3]] , 1+1,2+1]]    )         
  ], {\[Sigma],1,2},   {m,0,L1-1},{n,0,L2-1} ];     Chop[E/ Nc  ]      ];



(* ::Subsubsection::Bold::Closed:: *)
(*Mean field*)


(* ::Text::Bold:: *)
(*T matrix: from c Majorana fermions to f complex fermions *)


(* ::Code::Bold:: *)
(**)


Tx = KroneckerProduct[ {{1,1},{I,-I}}, SparseArray[{Band[{1,1}]-> {0,1,0,0}},{4,4}]  ]; 
Ty = KroneckerProduct[ {{1,1},{I,-I}}, SparseArray[{Band[{1,1}]-> {0,0,1,0}},{4,4}]  ];
Tz = KroneckerProduct[ {{1,1},{I,-I}}, SparseArray[{Band[{1,1}]-> {1,0,0,1}},{4,4}]  ]; 
Tmat[L1_,L2_] := Module[{Nc=L1 L2,T}, T=Sum[Module[{rx,ry,rz,mx,ny},
 mx=Mod[m+1,L1,1];ny=Mod[n+1,L2,1]; rx=mx+(n-1)L1;ry=m+(ny-1)L1;rz=m+(n-1)L1;
 KroneckerProduct[ one[rz,rx,Nc],Tx]+  KroneckerProduct[ one[rz,ry,Nc],Ty]+ KroneckerProduct[ one[rz,rz,Nc],Tz ]] ,{m,1,L1},{n,1,L2}];T];



(* ::Text::Bold:: *)
(*Eigen-values, eigen-vectors and U matrix*)


(* ::Code::Bold:: *)
(**)


Umat[H_]:= Module[ {R=Quiet@Eigensystem@N[H]},ReverseSort[R\[Transpose]]\[Transpose][[2]]\[Transpose] ];
EandU[H_]:= Module[ {R=Transpose@ReverseSort@Transpose@Quiet@Eigensystem@N[H]},{R[[1]],R[[2]]\[Transpose] }];
\[DoubleStruckCapitalU]occupied[TU_,Nc_]:= Drop[Take[TU\[Transpose],-4Nc-1],{2}]\[Transpose];


(* ::Text::Bold:: *)
(*Changing individual JK\[CapitalGamma] values due to electric field to trap the vortex*)


(* ::Code::Bold:: *)
(**)


uniform[K_,L1_,L2_]:=Table[   {K,K,K}     ,L1 L2]; (* for K = J,K,\[CapitalGamma]*)
addVortex[ Ko_,K0_, R_,L1_,L2_]:= Module[ {m=R[[1]],n=R[[2]],K=Ko,r=Table[0,6]},
r[[1]]=m+n L1 +1;r[[2]]=Mod[m,L1]+Mod[n+1,L2] L1 +1;r[[3]]=Mod[m-1,L1]+Mod[n+1,L2]L1 +1;r[[4]]=Mod[m-1,L1]+Mod[n+1,L2]L1 +1;r[[5]]=Mod[m-1,L1]+n L1 +1;r[[6]]=m+n L1 +1;
K[[ r[[1]],1]]=K0;K[[ r[[2]],3]]=K0;K[[ r[[3]],2]]=K0;
K[[ r[[4]],1]]=K0;K[[ r[[5]],3]]=K0;K[[ r[[6]],2]]=K0;
K];




(* ::Code::Bold:: *)
(**)


add2Vortices[ Ko_,K1_,K2_, R1_,R2_,L1_,L2_] :=Module[ {K}, K = addVortex[ Ko,K1, R1,L1,L2]; addVortex[ K,K2, R2,L1,L2]];
add2VorticesMaxSpaced[ Ko_,K1_,K2_,L1_,L2_] := Module[{R1,R2,d}, d=\[LeftFloor](L1-1)/2\[RightFloor]; R1= {\[LeftFloor]L1/2\[RightFloor] ,\[LeftFloor](L2+1)/2\[RightFloor] };R2= {\[LeftFloor]L1/2\[RightFloor]   +d,\[LeftFloor](L2+1)/2\[RightFloor] -d};
	add2Vortices[ Ko,K1,K2, R1,R2,L1,L2]];
flipBondsZMaxSpaced[\[Chi]0_,L1_,L2_]:= Module[{R1,d,\[Chi]=\[Chi]0,m0=\[LeftCeiling]L1/2\[RightCeiling] -1,n0=\[LeftCeiling](L2+1)/2\[RightCeiling]-1}, d=\[LeftFloor](L1-1)/2\[RightFloor]; R1= { m0,n0};
Do[  Module[{r,m,n}, m =Mod[R1[[1]]+i,L1] ; n=Mod[R1[[2]]-i,L2] ; r = m+n L1+1;       \[Chi][[3, r ]][[4,4]] =-\[Chi]0[[3, r ]][[4,4]] ;       \[Chi][[3, r ]][[1,1]] =-\[Chi]0[[3, r ]][[1,1]]  ]   , {i,0,d-1}];
\[Chi]           ];


(* ::Subsubsection::Bold::Closed:: *)
(*Mean field 2*)


(* ::Code::Bold:: *)
(**)


icc[U_,L_,T_]:=Module[  { Nc=L^2,TUh,icc}, TUh=Chop[(T . U)[[;;,-4Nc;;-1]],10^-12]; icc=I TUh . TUh\[ConjugateTranspose];icc=(icc-ConjugateTranspose[icc])/2   ];


(* ::Code::Bold:: *)
(**)


MFpParallel[U_,L_,T_] :=Module[  { Nc=L^2,TU,TUh,u,\[Chi],\[Omega]},
\[Chi]=Array[Null,{3,Nc,4,4}];\[Omega]=Array[Null,{2,Nc,4,4}];
TU=T . U;TUh=TU[[;;,-4Nc;;-1]];u=Chop[I TUh . TUh\[ConjugateTranspose],10^-12];
   Do[ Module[{m,n,\[Alpha],\[Beta],rz,rx,ry,Io}, rz=\[LeftFloor]R0/16\[RightFloor];\[Beta]=\[LeftFloor](R0-16rz)/4\[RightFloor]+1;\[Alpha]=R0-16rz-4(\[Beta]-1)+1;
   n=\[LeftFloor]rz/L\[RightFloor];m=rz-n L;rx=Mod[m+1,L]+n L;ry=m+Mod[n+1,L] L;Io=Mod[\[Alpha]+8rz,8Nc,1];
    \[Chi][[1,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8rx,8Nc,1],Io]];
    \[Chi][[2,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8ry,8Nc,1],Io]];
    \[Chi][[3,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8rz,8Nc,1],Io]];
    \[Omega][[1,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+8rz,8Nc,1],Io]];
    \[Omega][[2,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8rz,8Nc,1],Mod[Io+4,8Nc,1]] ];
    ]; , {R0,0,16Nc-1}  ];
{\[Chi][[1]],\[Chi][[2]],\[Chi][[3]],\[Omega][[1]],\[Omega][[2]],{{},{},{}},{{},{},{}},u }
];


(* ::Code::Bold:: *)
(**)


BiParallel[U_,L_,T_] :=Module[  { Nc=L^2,TU,TUh,\[Chi],\[Omega],\[Xi],u},
\[Chi]=Array[Null,{3,Nc,4,4}];\[Omega]=Array[Null,{2,Nc,4,4}];\[Xi]=Array[Null,{6,Nc,4,4}];
TU=T . U;TUh=TU[[;;,-4Nc;;-1]];u=Chop[I TUh . TUh\[ConjugateTranspose],10^-12];u=(u-ConjugateTranspose[u])/2;
Do[
Module[{rz,rx,ry,Io,R},rz=m+n L;rx=Mod[m+1,L]+n L;ry=m+Mod[n+1,L] L;Io=Mod[\[Alpha]+8rz,8Nc,1];
R={m+Mod[n+1,L]L,Mod[m-1,L]+n L,Mod[m+1,L]+Mod[n-1,L] L,m+Mod[n-1,L]L,Mod[m+1,L]+n L,Mod[m-1,L]+Mod[n+1,L] L};
\[Chi][[1,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8rx,8Nc,1],Mod[\[Alpha]+8rz,8Nc,1]]];
\[Chi][[2,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8ry,8Nc,1],Mod[\[Alpha]+8rz,8Nc,1]]];
\[Chi][[3,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8rz,8Nc,1],Mod[\[Alpha]+8rz,8Nc,1]]];
\[Omega][[1,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+8rz,8Nc,1],  Mod[\[Alpha]+8rz,8Nc,1]]];
\[Omega][[2,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8rz,8Nc,1],Mod[\[Alpha]+4+8rz,8Nc,1]]];

\[Xi][[1,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+8R[[1]],8Nc,1],Mod[\[Alpha]+8rz,8Nc,1]]];
\[Xi][[2,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+8R[[2]],8Nc,1],Mod[\[Alpha]+8rz,8Nc,1]]];
\[Xi][[3,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+8R[[3]],8Nc,1],Mod[\[Alpha]+8rz,8Nc,1]]];
\[Xi][[4,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8R[[4]],8Nc,1],Mod[\[Alpha]+4+8rz,8Nc,1]]];
\[Xi][[5,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8R[[5]],8Nc,1],Mod[\[Alpha]+4+8rz,8Nc,1]]];
\[Xi][[6,rz+1,\[Alpha],\[Beta]]]=u[[Mod[\[Beta]+4+8R[[6]],8Nc,1],Mod[\[Alpha]+4+8rz,8Nc,1]]];
]; 
, {n,0,L-1}, {m,0,L-1}, {\[Alpha],1,4},{\[Beta],1,4}  ];
{\[Chi][[1]],\[Chi][[2]],\[Chi][[3]],\[Omega][[1]],\[Omega][[2]],{\[Xi][[1]],\[Xi][[2]],\[Xi][[3]]},{\[Xi][[4]],\[Xi][[5]],\[Xi][[6]]} ,u }
];


(* ::Subsubsection::Bold::Closed:: *)
(*Mean field 3*)


(* ::Code::Bold:: *)
(**)


\[Chi]gauge4v[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;       r = m+n L+1;    \[Chi][[1,r]][[2,2]] =Abs@\[Chi]0[[1,r]][[2,2]]; \[Chi][[1,r]][[1,1]] =-Abs@\[Chi]0[[1,r]][[1,1]];      ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1; r = m+n L+1;    \[Chi][[1,r]][[2,2]] =Abs@\[Chi]0[[1,r]][[2,2]];   \[Chi][[1,r]][[1,1]] =-Abs@\[Chi]0[[1,r]][[1,1]];      ]   , {i,1,d1}];       \[Chi] ];

\[Chi]fixgauge[\[Chi]0_,u0_,L_]:= Module[{\[Chi]=\[Chi]0 }, 
Do[  \[Chi][[\[Alpha],r]][[\[Alpha]+1,\[Alpha]+1]] = u0[[r,\[Alpha]]]Abs@\[Chi]0[[\[Alpha],r]][[\[Alpha]+1,\[Alpha]+1]];      ,{\[Alpha],1,3} , {r,1,L^2}];   
   \[Chi] ];
\[Chi]gaugetransf[\[Chi]0_,u0_,L_]:= Module[{\[Chi]=\[Chi]0 }, 
Do[  \[Chi][[\[Alpha],r]][[1,1]] = -u0[[r,\[Alpha]]]\[Chi]0[[\[Alpha],r]][[1,1]]; \[Chi][[\[Alpha],r]][[\[Alpha]+1,\[Alpha]+1]] = -u0[[r,\[Alpha]]]\[Chi]0[[\[Alpha],r]][[\[Alpha]+1,\[Alpha]+1]];      ,{\[Alpha],1,3} , {r,1,L^2}];   
   \[Chi] ];
\[Chi]gaugetransf[\[Chi]0_,u0_,L_]:= Module[{\[Chi]=\[Chi]0 }, 
Do[  \[Chi][[\[Alpha],r]][[1,1]] = Abs@\[Chi]0[[\[Alpha],r]][[1,1]]; \[Chi][[\[Alpha],r]][[\[Alpha]+1,\[Alpha]+1]] = -Abs@\[Chi]0[[\[Alpha],r]][[\[Alpha]+1,\[Alpha]+1]];      ,{\[Alpha],1,3} , {r,1,L^2}];   
   \[Chi] ];
icc[U_,L_,T_]:=Module[  { Nc=L^2,TUh,icc}, TUh=Chop[(T . U)[[;;,-4Nc;;-1]],10^-12]; icc=I TUh . TUh\[ConjugateTranspose];(icc-ConjugateTranspose[icc])/2   ];

\[Chi]gauge4vCHANGE[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;       r = m+n L+1;    \[Chi][[1,r]][[2,2]] =-Abs@\[Chi]0[[1,r]][[2,2]];  \[Chi][[1,r]][[1,1]] =Abs@\[Chi]0[[1,r]][[1,1]];      ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1; r = m+n L+1;    \[Chi][[1,r]][[2,2]] =-Abs@\[Chi]0[[1,r]][[2,2]];\[Chi][[1,r]][[1,1]] =Abs@\[Chi]0[[1,r]][[1,1]];         ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =mS ; n=Mod[nS-i,L] ;       r = m+n L+1;    \[Chi][[1,r]][[2,2]] =Abs@\[Chi]0[[1,r]][[2,2]]; \[Chi][[1,r]][[1,1]] =-Abs@\[Chi]0[[1,r]][[1,1]];       ]   , {i,0,d1-1}];   
Do[  Module[{r,m,n}, m =mN ; n=Mod[nN+i+1,L]; r = m+n L+1;    \[Chi][[1,r]][[2,2]] =Abs@\[Chi]0[[1,r]][[2,2]];   \[Chi][[1,r]][[1,1]] =-Abs@\[Chi]0[[1,r]][[1,1]];      ]   , {i,0,d1-1}];       \[Chi] ];

TranslationMF2[\[Chi]0_,\[Omega]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,\[Omega]=\[Omega]0,mS,nS,mN,nN,r2,Nc=L^2},   
r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[ Module[{m,n,mT,nT}, 
   n=\[LeftFloor]r/L\[RightFloor];m=r-n L; nT=Mod[n+2r2,L];
    \[Chi][[1,m+L n+1]]=\[Chi]0[[1]][[m+L nT+1]] ;
    \[Chi][[2,m+L n+1]]=\[Chi]0[[2]][[m+L nT+1]] ;
    \[Chi][[3,m+L n+1]]=\[Chi]0[[3]][[m+L nT+1]] ;
    \[Omega][[1,m+L n+1]]=\[Omega]0[[1]][[m+L nT+1]] ;
    \[Omega][[2,m+L n+1]]=\[Omega]0[[2]][[m+L nT+1]] ;
]; , {r,0,Nc-1}  ];
{\[Chi],\[Omega]}];


(* ::Code::Bold:: *)
(**)


(*\[Chi]gauge4vChangeXtoY[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;          r=m+n L+1;   \[Chi][[1,r]][[2,2]] =-Abs@\[Chi]0[[1,r]][[2,2]];   \[Chi][[1,r]][[1,1]] =Abs@\[Chi]0[[1,r]][[1,1]];      ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1;        r=m+n L+1;    \[Chi][[1,r]][[2,2]] =-Abs@\[Chi]0[[1,r]][[2,2]];  \[Chi][[1,r]][[1,1]] =Abs@\[Chi]0[[1,r]][[1,1]];         ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =mS+i ; n=nS ;         r=m+n L+1;    \[Chi][[2,r]][[3,3]] =Abs@\[Chi]0[[2,r]][[3,3]];   \[Chi][[2,r]][[1,1]] =-Abs@\[Chi]0[[2,r]][[1,1]];       ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =mN -i+1 ; n=nN;        r=m+n L+1;    \[Chi][[2,r]][[3,3]] =Abs@\[Chi]0[[2,r]][[3,3]];   \[Chi][[2,r]][[1,1]] =-Abs@\[Chi]0[[2,r]][[1,1]];      ]   , {i,1,d1}];       \[Chi] ];

\[Chi]gauge4vChangeXtoZ[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,mW,nW,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;mW=r2-1;nW=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;          r=m+n L+1;   \[Chi][[1,r]][[2,2]] =-Abs@\[Chi]0[[1,r]][[2,2]];   \[Chi][[1,r]][[1,1]] =Abs@\[Chi]0[[1,r]][[1,1]];      ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1;        r=m+n L+1;    \[Chi][[1,r]][[2,2]] =-Abs@\[Chi]0[[1,r]][[2,2]];  \[Chi][[1,r]][[1,1]] =Abs@\[Chi]0[[1,r]][[1,1]];         ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =Mod[mS+i,L]; n=Mod[nS-i+1,L]; r=m+n L+1;    \[Chi][[3,r]][[4,4]] =Abs@\[Chi]0[[3,r]][[4,4]];   \[Chi][[3,r]][[1,1]] =-Abs@\[Chi]0[[3,r]][[1,1]];       ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =Mod[mW+i,L]; n=Mod[nW-i+1,L]; r=m+n L+1;    \[Chi][[3,r]][[4,4]] =Abs@\[Chi]0[[3,r]][[4,4]];   \[Chi][[3,r]][[1,1]] =-Abs@\[Chi]0[[3,r]][[1,1]];       ]   , {i,1,d1}];       \[Chi] ]; *)


(* ::Code::Bold:: *)
(**)


(*changeSign[m_,s_]:={{s m\[LeftDoubleBracket]1,1\[RightDoubleBracket],m\[LeftDoubleBracket]1,2\[RightDoubleBracket],m\[LeftDoubleBracket]1,3\[RightDoubleBracket],m\[LeftDoubleBracket]1,4\[RightDoubleBracket]},{m\[LeftDoubleBracket]2,1\[RightDoubleBracket],s m\[LeftDoubleBracket]2,2\[RightDoubleBracket],s m\[LeftDoubleBracket]2,3\[RightDoubleBracket],s m\[LeftDoubleBracket]2,4\[RightDoubleBracket]},{m\[LeftDoubleBracket]3,1\[RightDoubleBracket],s m\[LeftDoubleBracket]3,2\[RightDoubleBracket],s m\[LeftDoubleBracket]3,3\[RightDoubleBracket],s m\[LeftDoubleBracket]3,4\[RightDoubleBracket]},{m\[LeftDoubleBracket]4,1\[RightDoubleBracket],s m\[LeftDoubleBracket]4,2\[RightDoubleBracket],s m\[LeftDoubleBracket]4,3\[RightDoubleBracket],s m\[LeftDoubleBracket]4,4\[RightDoubleBracket]}};
\[Chi]gauge4vChangeX[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;         r=m+n L+1;   \[Chi][[1,r]]=changeSign[\[Chi]0[[1,r]],-1];  \[Chi][[2,r]]=-changeSign[\[Chi]0[[2,r]],-1]; \[Chi][[3,r]]=-changeSign[\[Chi]0[[3,r]],-1];         ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1;        r=m+n L+1;   \[Chi][[1,r]]=changeSign[\[Chi]0[[1,r]],-1];  \[Chi][[2,r]]=-changeSign[\[Chi]0[[2,r]],-1]; \[Chi][[3,r]]=-changeSign[\[Chi]0[[3,r]],-1];           ]   , {i,1,d1}]; \[Chi] ];
\[Chi]gauge4vChangeXtoY[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;         r=m+n L+1;   \[Chi][[1,r]] =changeSign[\[Chi]0[[1,r]],-1];          ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1;        r=m+n L+1;   \[Chi][[1,r]] =changeSign[\[Chi]0[[1,r]],-1];          ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =mS+i ; n=nS ;         r=m+n L+1;   \[Chi][[2,r]] =changeSign[\[Chi]0[[2,r]],-1];          ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =mN -i+1 ; n=nN;       r=m+n L+1;   \[Chi][[2,r]] =changeSign[\[Chi]0[[2,r]],-1];          ]   , {i,1,d1}];       \[Chi] ];

\[Chi]gauge4vChangeXtoZ[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,mW,nW,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;mW=r2-1;nW=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;                     r=m+n L+1;     \[Chi][[1,r]]=changeSign[\[Chi]0[[1,r]],-1];      ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1;                    r=m+n L+1;     \[Chi][[1,r]]=changeSign[\[Chi]0[[1,r]],-1];          ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =Mod[mS+i,L]; n=Mod[nS-i+1,L];   r=m+n L+1;     \[Chi][[3,r]]=changeSign[\[Chi]0[[3,r]],-1];        ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =Mod[mW+i,L]; n=Mod[nW-i+1,L];   r=m+n L+1;     \[Chi][[3,r]]=changeSign[\[Chi]0[[3,r]],-1];        ]   , {i,1,d1}];       \[Chi] ];*)


(* ::Code::Bold:: *)
(**)


changeSign[m_,s_]:={{s m[[1,1]],m[[1,2]],m[[1,3]],m[[1,4]]},{m[[2,1]],s m[[2,2]],s m[[2,3]],s m[[2,4]]},{m[[3,1]],s m[[3,2]],s m[[3,3]],s m[[3,4]]},{m[[4,1]],s m[[4,2]],s m[[4,3]],s m[[4,4]]}};
\[Chi]gauge4vChangeX[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;         r=m+n L+1;   \[Chi][[1,r]]=-\[Chi]0[[1,r]];          ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1;        r=m+n L+1;   \[Chi][[1,r]]=-\[Chi]0[[1,r]];           ]   , {i,1,d1}]; \[Chi] ];
\[Chi]gauge4vChangeXtoY[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;         r=m+n L+1;   \[Chi][[1,r]]=-\[Chi]0[[1,r]];           ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1;        r=m+n L+1;   \[Chi][[1,r]]=-\[Chi]0[[1,r]];           ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =mS+i ; n=nS ;         r=m+n L+1;   \[Chi][[2,r]]=-\[Chi]0[[2,r]];           ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =mN -i+1 ; n=nN;       r=m+n L+1;   \[Chi][[2,r]]=-\[Chi]0[[2,r]];           ]   , {i,1,d1}];       \[Chi] ];

\[Chi]gauge4vChangeXtoZ[\[Chi]0_,L_]:= Module[{d1,d2,\[Chi]=\[Chi]0,mS,nS,mN,nN,mW,nW,r2 },   r2= \[LeftFloor]1/2 \[LeftCeiling]L/2\[RightCeiling]\[RightFloor]; d2=2r2;d1=L-d2;  mS=r2-1;nS=r2-1;mN=L-r2-1;nN=L-r2-1;mW=r2-1;nW=L-r2-1;
Do[  Module[{r,m,n}, m =mS ; n=nS+i ;                     r=m+n L+1;     \[Chi][[1,r]]=-\[Chi]0[[1,r]];       ]   , {i,1,d1}];   
Do[  Module[{r,m,n}, m =mN ; n=nN-i+1;                    r=m+n L+1;     \[Chi][[1,r]]=-\[Chi]0[[1,r]];           ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =Mod[mS+i,L]; n=Mod[nS-i+1,L];   r=m+n L+1;     \[Chi][[3,r]]=-\[Chi]0[[3,r]];         ]   , {i,1,d1}];
Do[  Module[{r,m,n}, m =Mod[mW+i,L]; n=Mod[nW-i+1,L];   r=m+n L+1;     \[Chi][[3,r]]=-\[Chi]0[[3,r]];        ]   , {i,1,d1}];       \[Chi] ];


(* ::Code::Bold:: *)
(**)


reflectXinternal[m_]:={
{m[[1,1]],m[[1,3]],m[[1,2]],m[[1,4]]},
{m[[3,1]],m[[3,3]],m[[3,2]],m[[3,4]]},
{m[[2,1]],m[[2,3]],m[[2,2]],m[[2,4]]},
{m[[4,1]],m[[4,3]],m[[4,2]],m[[4,4]]} }; 
reflectX[\[Chi]0_,\[Omega]0_,L_]:=Module[{\[Chi]=\[Chi]0,\[Omega]=\[Omega]0,Nc=L^2},
Do[ Module[{m,n}, 
   n=\[LeftFloor]r/L\[RightFloor];m=r-n L;(*rx=Mod[m+1,L]+n L;ry=m+Mod[n+1,L] L;Io=Mod[\[Alpha]+8rz,8Nc,1];rz=\[LeftFloor]R0/16\[RightFloor];\[Beta]=\[LeftFloor]1+(R0-16rz)/4\[RightFloor];\[Alpha]=R0-16rz-4(\[Beta]-1)+1; *)
    \[Chi][[1,m+L n +1]]=reflectXinternal@\[Chi]0[[2]][[n+L m +1]] ;
    \[Chi][[2,m+L n +1]]=reflectXinternal@\[Chi]0[[1]][[n+L m +1]] ;
    \[Chi][[3,m+L n +1]]=reflectXinternal@\[Chi]0[[3]][[n+L m +1]] ;
    \[Omega][[1,m+L n +1]]=reflectXinternal@\[Omega]0[[1]][[n+L m +1]] ;
    \[Omega][[2,m+L n +1]]=reflectXinternal@\[Omega]0[[2]][[n+L m +1]] ;
    ]; , {r,0,Nc-1}  ];   
{\[Chi],\[Omega]}];


(* ::Subsubsection::Bold::Closed:: *)
(*Couplings  *)


(* ::Code::Bold:: *)
(**)


Jc[eV0_,JH_,U_,t_]:=1/54 ((2 (t[[1]]-t[[3]])^2)/(eV0-JH+U)-(2 (t[[1]]-t[[3]])^2)/(eV0+JH-U)-(6 t[[1]] (t[[1]]+2 t[[3]]))/(eV0+3 JH-U)+(6 t[[1]] (t[[1]]+2 t[[3]]))/(eV0-3 JH+U)+(2 t[[1]]+t[[3]])^2/(-eV0+2 JH+U)+(2 t[[1]]+t[[3]])^2/(eV0+2 JH+U));
Kc[eV0_,JH_,U_,t_]:=(2JH )/9 ( (t[[1]]-t[[3]])^2-3 t[[2]]^2 ) (eV0^2+3 JH^2-4 JH U+U^2)/((eV0+JH-U) (eV0+3 JH-U) (eV0-3 JH+U) (eV0-JH+U));
\[CapitalGamma]c[eV0_,JH_,U_,t_]:= (4 JH t[[2]] (t[[1]]-t[[3]]) )/9  (eV0^2+3 JH^2-4 JH U+U^2)/( (eV0+JH-U) (eV0+3 JH-U) (eV0-3 JH+U) (eV0-JH+U));


(* ::Code::Bold:: *)
(**)


Jr[eV0_,JH_,U_,t_]:=Jc[eV0,JH,U,t]/Abs[Kc[0,JH,U,t]];
Kr[eV0_,JH_,U_,t_]:=Kc[eV0,JH,U,t]/Abs[Kc[0,JH,U,t]];
\[CapitalGamma]r[eV0_,JH_,U_,t_]:= \[CapitalGamma]c[eV0,JH,U,t]/Abs[Kc[0,JH,U,t]];


(* ::Subsection::Bold::Closed:: *)
(*Momentum def*)


(* ::Subsubsection::Bold::Closed:: *)
(*momentum*)


(* ::Code::Bold:: *)
(**)


mx= 2\[Pi] {1 , 1/Sqrt[3]}; my= 2\[Pi] {-1 , 1/Sqrt[3]};


(* ::Code::Bold:: *)
(**)


(*toMomentum[n_,Nb_]:= If[FractionalPart[Sqrt@Nb] \[NotEqual] 0, Print["Nb=",Nb,", is not a perfect square."];Abort[],
				Module[ {L = IntegerPart[Sqrt@Nb],p1,p2},  p1 =  Mod[n-1,L]+1; p2 =   (n-p1)/L+1;   mxp1/L +myp2/L]
];*)
toMomentum[n_,Nb_]:= Module[ {L = IntegerPart[Sqrt@Nb],p1,p2},  p1 =  Mod[n,L,1]; p2 =   (n-p1)/L+1;   mx p1/L +my p2/L];
toMomentumTable[L_]:=Table[Module[{p1,p2},p1=Mod[n,L,1];p2=(n-p1)/L+1;{(2(p1-p2)\[Pi])/L,(2 (p1+p2) \[Pi])/(Sqrt[3] L)}],{n,1,L^2}];
toMomentumInverse[k_,Nb_]:= If[FractionalPart[Sqrt@Nb] != 0, Print["Nb=",Nb,", is not a perfect square."];Abort[],    
Module[ {M = IntegerPart[Sqrt@Nb],p1,p2},   p1=k . nx M/(2\[Pi]); p2=k . ny M/(2\[Pi]);  Round[p1 + (p2-1)M ]  ]];
kD2=1/3 mx+2/3 my;
kD1=2/3 mx+1/3 my;


(* ::Subsubsection::Bold::Closed:: *)
(*HMF*)


(* ::Code::Bold:: *)
(**)


HMFk[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_,\[Eta]_, k_] := 
Module[{hx,hy,hz,kx,ky,H},
kx=k . nx;hx=Re@Hx[J,K,\[CapitalGamma],h,\[Chi],\[Omega]]; hy=Re@Hy[J,K,\[CapitalGamma],h,\[Chi],\[Omega]];hz=Re@Hz[J,K,\[CapitalGamma],h,\[Chi],\[Omega],\[Eta]];ky=k . ny; H=I(hx) Exp[I kx]+I (hy)Exp[I ky]+I(hz); N@(H+H\[ConjugateTranspose] ) ];

Tk=KroneckerProduct[ {{1,1},{I,-I}}, SparseArray[{Band[{1,1}]-> {1,1,1,1}},{4,4}] ]; 
UmatK[H_]:= Module[ {R=Quiet@Eigensystem@N[H]},ReverseSort[R\[Transpose]]\[Transpose][[2]]\[Transpose] ];

MFparK[J_,K_,\[CapitalGamma]_,h_,\[Eta]_, \[Chi]0_,\[Omega]0_,L_]:=Module[ { Nc=L^2,u,ux,uy,uz,\[Chi]=Array[Null,3],\[Omega]=Array[Null,2]},
u=Total@Table[ Module[{H,Hr,U,TU,TUh,k,uu},
k=toMomentum[l,Nc];
H=HMFk[J,K,\[CapitalGamma],h,\[Chi]0,\[Omega]0,\[Eta], k];Hr=Tk\[ConjugateTranspose] . H . Tk; U=UmatK[Hr];
TU=Tk . U;TUh=TU[[;;,-4;;-1]];uu=(1/Nc)Chop[I TUh . TUh\[ConjugateTranspose] ];uu=(uu-Transpose[uu])/2;
{uu Exp[I k . nx],uu Exp[I k . ny],uu }  ],{l,1,Nc}]; 
\[Chi][[1]]=Table[ u[[1]][[\[Beta]+4,\[Alpha]]]   ,  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Chi][[2]]=Table[ u[[2]][[\[Beta]+4,\[Alpha]]]   ,  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Chi][[3]]=Table[ u[[3]][[\[Beta]+4,\[Alpha]]]   ,  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Omega][[1]]=Table[ u[[3]][[\[Beta],\[Alpha]]]     ,  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Omega][[2]]=Table[ u[[3]][[\[Beta]+4,\[Alpha]+4]] ,  {\[Alpha],1,4},{\[Beta],1,4} ];
{\[Chi][[1]],\[Chi][[2]],\[Chi][[3]],\[Omega][[1]],\[Omega][[2]]}
];



(* ::Subsubsection::Bold::Closed:: *)
(*Auxiliary matrices for the Hamiltonian  [8x8 matrices]*)


(* ::Code::Bold:: *)
(**)


Hx[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_]:={
{0,0,0,0,-K \[Chi][[1,2,2]],K \[Chi][[1, 2,1]],0,0},
{0,0,0,0, K \[Chi][[1, 1,2]],-K \[Chi][[1, 1,1]],0,0},
{0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,J (-\[Chi][[1, 2,2]]-\[Chi][[1, 3,3]]-\[Chi][[1, 4,4]]),J \[Chi][[1, 2,1]],J \[Chi][[1, 3,1]],J \[Chi][[1, 4,1]]},
{0,0,0,0,J \[Chi][[1, 1,2]],-J \[Chi][[1, 1,1]],0,0},
{0,0,0,0,J \[Chi][[1, 1,3]],0,-J \[Chi][[1, 1,1]],0},
{0,0,0,0,J \[Chi][[1, 1,4]],0,0,-J \[Chi][[1, 1,1]]},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,\[CapitalGamma] (-\[Chi][[1,3,4]]-\[Chi][[1, 4,3]]),0,\[CapitalGamma] \[Chi][[1, 4,1]],\[CapitalGamma] \[Chi][[1, 3,1]]},
{0,0,0,0,0,0,0,0},
{0,0,0,0,\[CapitalGamma] \[Chi][[1, 1,4]],0,0,\[CapitalGamma] \[Chi][[1, 1,1]]},
{0,0,0,0,\[CapitalGamma] \[Chi][[1, 1,3]],0,\[CapitalGamma] \[Chi][[1, 1,1]],0},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}; 
Hy[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_]:={
{0,0,0,0,-K \[Chi][[2, 3,3]],0,K \[Chi][[2, 3,1]],0},
{0,0,0,0,0,0,0,0},
{0,0,0,0,K \[Chi][[2, 1,3]],0,-K \[Chi][[2, 1,1]],0},
{0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,J (-\[Chi][[2, 2,2]]-\[Chi][[2, 3,3]]-\[Chi][[2, 4,4]]),J \[Chi][[2, 2,1]],J \[Chi][[2, 3,1]],J \[Chi][[2, 4,1]]},
{0,0,0,0,J \[Chi][[2, 1,2]],-J \[Chi][[2, 1,1]],0,0},
{0,0,0,0,J \[Chi][[2, 1,3]],0,-J \[Chi][[2, 1,1]],0},
{0,0,0,0,J \[Chi][[2, 1,4]],0,0,-J \[Chi][[2, 1,1]]},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,\[CapitalGamma] (-\[Chi][[2, 2,4]]-\[Chi][[2, 4,2]]),\[CapitalGamma] \[Chi][[2, 4,1]],0,\[CapitalGamma] \[Chi][[2, 2,1]]},
{0,0,0,0,\[CapitalGamma] \[Chi][[2, 1,4]],0,0,\[CapitalGamma] \[Chi][[2, 1,1]]},
{0,0,0,0,0,0,0,0},
{0,0,0,0,\[CapitalGamma] \[Chi][[2, 1,2]],\[CapitalGamma] \[Chi][[2, 1,1]],0,0},
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}};

Hz[J_,K_,\[CapitalGamma]_,h_,\[Chi]_,\[Omega]_,\[Eta]_]:=Module[{\[Lambda]1,\[Lambda]2,m,\[Delta]\[Omega]},
\[Delta]\[Omega]=Table[{ \[Omega][[\[Sigma],2,1]]-\[Omega][[\[Sigma],3,4]],\[Omega][[\[Sigma],3,1]]-\[Omega][[\[Sigma],4,2]],\[Omega][[\[Sigma],4,1]]-\[Omega][[\[Sigma],2,3]]}+I 10^(-9.9)  ,{\[Sigma],1,2}];
m=Table[ \[Omega][[\[Sigma],2,1]]\[Omega][[\[Sigma],3,4]]+\[Omega][[\[Sigma],3,1]]\[Omega][[\[Sigma],4,2]]+\[Omega][[\[Sigma],4,1]]\[Omega][[\[Sigma],2,3]] +1 ,{\[Sigma],1,2}];
\[Lambda]2=Re@Table[{\[Omega][[\[Sigma],2,1]]/\[Delta]\[Omega][[\[Sigma],1]],\[Omega][[\[Sigma],3,1]]/\[Delta]\[Omega][[\[Sigma],2]] ,\[Omega][[\[Sigma],4,1]]/\[Delta]\[Omega][[\[Sigma],3]]  },{\[Sigma],1,2}];
\[Lambda]1=Re@Table[{-\[Omega][[\[Sigma],3,4]]/\[Delta]\[Omega][[\[Sigma],1]] ,-\[Omega][[\[Sigma],4,2]]/\[Delta]\[Omega][[\[Sigma],2]] ,-\[Omega][[\[Sigma],2,3]]/\[Delta]\[Omega][[\[Sigma],3]] },{\[Sigma],1,2}];
{
{0,0,0,0,         -K \[Chi][[3, 4,4]],0,0, K \[Chi][[3, 4,1]]    },
{0,0,0,0,         0,0,0,0   },
{0,0,0,0,         0,0,0,0   },
{0,0,0,0,         K \[Chi][[3, 1,4]],0,0,-K \[Chi][[3, 1,1]]  },
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,         J(-\[Chi][[3, 2,2]]-\[Chi][[3, 3,3]]-\[Chi][[3, 4,4]]),J \[Chi][[3, 2,1]],J \[Chi][[3, 3,1]],J \[Chi][[3, 4,1]]    },
{0,0,0,0,         J \[Chi][[3, 1,2]],-J \[Chi][[3, 1,1]],0,0   },
{0,0,0,0,         J \[Chi][[3, 1,3]],0,-J \[Chi][[3, 1,1]], 0  },
{0,0,0,0,         J \[Chi][[3, 1,4]],0,0,-J \[Chi][[3, 1,1]]  },
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,         \[CapitalGamma](-\[Chi][[3, 2,3]]-\[Chi][[3, 3,2]]),\[CapitalGamma] \[Chi][[3, 3,1]],\[CapitalGamma] \[Chi][[3, 2,1]],0   },
{0,0,0,0,         \[CapitalGamma] \[Chi][[3, 1,3]],0,\[CapitalGamma] \[Chi][[3, 1,1]],0   },
{0,0,0,0,         \[CapitalGamma] \[Chi][[3, 1,2]],\[CapitalGamma] \[Chi][[3, 1,1]],0, 0  },
{0,0,0,0,         0,0,0,0  },
{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}+{
{0,0,0,0,             0,0,0,0 },
{-heffM0[J,K,\[CapitalGamma],h,\[Omega], 2][[1]] + \[Lambda]1[[2,1]] heffM[J,K,\[CapitalGamma],h,\[Omega], 2][[1]],0, \[Lambda]2[[2,3]] heffM[J,K,\[CapitalGamma],h,\[Omega], 2][[3]],0,         0,0,0,0   },
{-heffM0[J,K,\[CapitalGamma],h,\[Omega], 2][[2]] + \[Lambda]1[[2,2]] heffM[J,K,\[CapitalGamma],h,\[Omega], 2][[2]],0,0,\[Lambda]2[[2,1]]  heffM[J,K,\[CapitalGamma],h,\[Omega], 2][[1]],         0,0,0,0   },
{-heffM0[J,K,\[CapitalGamma],h,\[Omega], 2][[3]] + \[Lambda]1[[2,3]]heffM[J,K,\[CapitalGamma],h,\[Omega], 2][[3]],\[Lambda]2[[2,2]] heffM[J,K,\[CapitalGamma],h,\[Omega], 2][[2]],0,0,         0,0,0,0   },
{0,0,0,0,0,0,0,0},
{0,0,0,0,-heffM0[J,K,\[CapitalGamma],h,\[Omega], 1][[1]] + \[Lambda]1[[1,1]] heffM[J,K,\[CapitalGamma],h,\[Omega], 1][[1]],0, \[Lambda]2[[1,3]] heffM[J,K,\[CapitalGamma],h,\[Omega], 1][[3]],0},
{0,0,0,0,-heffM0[J,K,\[CapitalGamma],h,\[Omega], 1][[2]] + \[Lambda]1[[1,2]]heffM[J,K,\[CapitalGamma],h,\[Omega], 1][[2]],0,0, \[Lambda]2[[1,1]] heffM[J,K,\[CapitalGamma],h,\[Omega], 1][[1]] },
{0,0,0,0,-heffM0[J,K,\[CapitalGamma],h,\[Omega], 1][[3]] + \[Lambda]1[[1,3]]heffM[J,K,\[CapitalGamma],h,\[Omega], 1][[3]], \[Lambda]2[[1,2]]heffM[J,K,\[CapitalGamma],h,\[Omega], 1][[2]],0,0} }];


(* ::Subsubsection::Bold::Closed:: *)
(*heff momentum*)


(* ::Code::Bold:: *)
(**)


heffM0[J_,K_,\[CapitalGamma]_,h_,\[Omega]_,\[Sigma]_]:={ 
h[[1]] -  K \[Omega][[\[Sigma], 2,1]]-3J \[Omega][[\[Sigma], 2,1]] -2\[CapitalGamma] (\[Omega][[\[Sigma], 3,1]] +\[Omega][[\[Sigma], 4,1]]), 
h[[2]] -  K \[Omega][[\[Sigma], 3,1]]-3 J \[Omega][[\[Sigma], 3,1]]-2\[CapitalGamma] (\[Omega][[\[Sigma], 2,1]] +\[Omega][[\[Sigma], 4,1]]), 
h[[3]] -  K \[Omega][[\[Sigma], 4,1]]-3 J \[Omega][[\[Sigma], 4,1]]-2\[CapitalGamma] (\[Omega][[\[Sigma], 2,1]] +\[Omega][[\[Sigma], 3,1]])};
heffM[J_,K_,\[CapitalGamma]_,h_,\[Omega]_,\[Sigma]_]:=heffM0[J,K,\[CapitalGamma],h,\[Omega],\[Sigma]];c
sym\[Omega][\[Omega]_]:= Module[ {\[Omega]A,\[Omega]B},
\[Omega]A=(1/2){  {I,0,0,0},  
{\[Omega][[1]][[2,1]]-\[Omega][[1]][[3,4]],I,-\[Omega][[1]][[4,1]]+\[Omega][[1]][[2,3]],0},  
{\[Omega][[1]][[3,1]]-\[Omega][[1]][[4,2]],0,I,-\[Omega][[1]][[2,1]]+\[Omega][[1]][[3,4]]},  
{\[Omega][[1]][[4,1]]-\[Omega][[1]][[2,3]],-\[Omega][[1]][[3,1]]+\[Omega][[1]][[4,2]],0,I}  };
\[Omega]B=(1/2){  {I,0,0,0},  
{\[Omega][[2]][[2,1]]-\[Omega][[2]][[3,4]],I,-\[Omega][[2]][[4,1]]+\[Omega][[2]][[2,3]],0},  
{\[Omega][[2]][[3,1]]-\[Omega][[2]][[4,2]],0,I,-\[Omega][[2]][[2,1]]+\[Omega][[2]][[3,4]]},  
{\[Omega][[2]][[4,1]]-\[Omega][[2]][[2,3]],-\[Omega][[2]][[3,1]]+\[Omega][[2]][[4,2]],0,I}  };
{\[Omega]A-ConjugateTranspose[\[Omega]A],\[Omega]B-ConjugateTranspose[\[Omega]B]}  ];
\[Omega]GA = {{I,0.000001,0.000001,0.000001},{-0.000001,I,0.000001,-0.000001},{-0.000001,-0.000001,I,0.000001},{-0.000001,0.000001,-0.000001,I}}; \[Omega]GB = \[Omega]GA;


(* ::Section::Bold:: *)
(*MF Loop -- microscopic parameters*)


ts = Table[ {5x,160,-12x,0,-60},{x,tV}];
hs =Table[  h[[1]]  hAngle[h[[2]],h[[3]]] , {h,hV}];  
eV0=0;U=2600;JH=300;
(*  \[Kappa]=0.1 : h=0.3292;   \[Kappa]=0.2 : 0.41475; \[Kappa]=0.05 : 0.2612 *)
parameters=Table[Flatten[ Table[ {N@Jr[0,JH,U,ts[[t]] ],N@Kr[0,JH,U,ts[[t]]],N@\[CapitalGamma]r[0,JH,U,ts[[t]]],hs[[h]] ,N@Jr[eVs[[ev]],JH,U,ts[[t]] ],N@Kr[eVs[[ev]],JH,U,ts[[t]]],N@\[CapitalGamma]r[eVs[[ev]],JH,U,ts[[t]]],hV[[h]],tV[[t]],eVs[[ev]]} , {t,1,Length@tV},  {h,1,Length@hV}],1] ,  {ev,1,Length@eVs}];


Print[" "];
Print["    NbName=",NbName,"; "];
Print["    Ls=",Ls,"; "];
Print["    tV=",tV,";    ts=", ts,"; "];
Print["    hV=",hV,";"]
Print["    hs=", hs,"; "];
Print["    Steps=",steps,"; "];
Print["    acuracy=",acuracy,"; "];
Print["    eVs=",eVs,"; "];
Print[" "];
Print["    Parameters="];
Do[
Print[ parameters[[i,j]] ],
{i,1,Length@parameters},
{j,1,Length@parameters[[i]] }
];
Print[" "];


(* ::Subsubsection::Bold:: *)
(*vortex free*)


tvf=AbsoluteTime[];Print["Definition timing= ",round[tvf-t0] ]; t0=tvf;


Print[" "];Print[" "];Print["    Starting free loop"];Print[" "];
t0=AbsoluteTime[]; 
Do[ Do[ Module[{ loaddata,\[CapitalGamma],J,K,L=Ls[[l]],Nc,h,\[CapitalLambda],T,H,\[Xi],EnG0,En,EnList={{},{},{}},u,u2,\[Chi],\[Omega],j,\[CapitalDelta]1=1,\[CapitalDelta]2=1,ES,gap,\[CapitalDelta]t,kTable,\[CapitalDelta]\[Omega],\[CapitalDelta]\[Omega]seq={},\[CapitalDelta]seq={},EMF,Esum,E\[Lambda],\[Eta]=\[Lambda]0,hp=Mod[p,Length@hV,1]},
{J,K,\[CapitalGamma],h}=parameters[[1,p]][[1;;4]]; Nc=L^2;
If[ p==1, \[Chi]G={\[Chi]Gx,\[Chi]Gy,\[Chi]Gz}; \[Omega]G={\[Omega]GA,\[Omega]GB}; ];\[Chi]=\[Chi]G; \[Omega]=\[Omega]G;   (*Print[MatrixForm/@\[Chi],MatrixForm/@\[Omega]];*)
Print["J=",J,"; K=",K, "; G=",\[CapitalGamma], "; L=",L,"; h=(", hV[[ hp,1 ]],",",hV[[ hp,2 ]],",",hV[[ hp,3]],"); "];
kTable=toMomentumTable[L];
For[j=1,( ( j<steps)\[And](Chop[ \[CapitalDelta]1,10^(-acuracy-6) ]!= 0) ), j++,  
If[j<=1, loaddata=loadDataTry[toPath[parameters[[1,p]],L,acuracy,"free",NbName]  ];
If[ !(loaddata===$Failed),{j,L,\[Chi],\[Omega],\[Xi],EnG0}=loaddata]];
u=Chop@Total@Table[ Module[{H0,Hr,U,TU,TUh,k,uu},
k=kTable[[l]];
H0=HMFk[J,K,\[CapitalGamma],h,\[Chi],\[Omega],\[Eta],k];
Hr=Tk\[ConjugateTranspose] . H0 . Tk; 
U=UmatK[Hr];
TU=Tk . U;
TUh=TU[[;;,-4;;-1]];
uu=(1/Nc)I TUh . TUh\[ConjugateTranspose] ;
uu=(uu-ConjugateTranspose[uu])/2;  {uu Exp[I k . nx],uu Exp[I k . ny],uu }  ],{l,1,Nc}]; 
\[Chi][[1]]=Table[ u[[1]][[\[Beta]+4,\[Alpha]]]   ,  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Chi][[2]]=Table[ u[[2]][[\[Beta]+4,\[Alpha]]]   ,  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Chi][[3]]=Table[ u[[3]][[\[Beta]+4,\[Alpha]]]   ,  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Omega][[1]]=Table[ u[[3]][[\[Beta],\[Alpha]]]     ,  {\[Alpha],1,4},{\[Beta],1,4} ];
\[Omega][[2]]=Table[ u[[3]][[\[Beta]+4,\[Alpha]+4]] ,  {\[Alpha],1,4},{\[Beta],1,4} ];
If[j>=2,\[CapitalDelta]2=\[CapitalDelta]1; \[CapitalDelta]1=Max[ Abs@(u-u2) ];\[CapitalDelta]seq={\[CapitalDelta]seq,{j,\[CapitalDelta]1}}  ;  
\[CapitalDelta]\[Omega]=1/6 Sum[Abs[\[Omega][[\[Sigma],2,1]]+\[Omega][[\[Sigma],3,4]]]+Abs[\[Omega][[\[Sigma],3,1]]+\[Omega][[\[Sigma],4,2]]]+Abs[\[Omega][[\[Sigma],4,1]]+\[Omega][[\[Sigma],2,3]]],{\[Sigma],1,2}]; \[CapitalDelta]\[Omega]seq={\[CapitalDelta]\[Omega]seq,{j,\[CapitalDelta]\[Omega]}} ;]; 
u2=u;  (*
If[j\[LessEqual]3\[Or]j\[GreaterEqual]steps-2,Print["j=",j,"; \[CapitalDelta]=",\[CapitalDelta]1,"; \[CapitalDelta]\[Omega]=",\[CapitalDelta]\[Omega], ";  ",MatrixForm/@{  \[Chi][[1]],\[Omega][[1]] }];  ] *)
	];
t1=AbsoluteTime[]; \[CapitalDelta]t=UnitConvert[ Quantity[t1 -t0, "Seconds" ], "Minutes" ]; t0=t1;
	\[CapitalDelta]seq=Partition[Flatten[\[CapitalDelta]seq],2]; \[CapitalDelta]\[Omega]seq=Partition[Flatten[\[CapitalDelta]\[Omega]seq],2];
Module[{\[Omega]v,\[Chi]v},
\[Omega]v={Table[\[Omega][[1]],{r,1,Nc} ],Table[\[Omega][[2]],{r,1,Nc} ] };
\[Chi]v={Table[\[Chi][[1]],{r,1,Nc} ], Table[\[Chi][[2]],{r,1,Nc}], Table[\[Chi][[3]],{r,1,Nc}]  };
EMF=EnMF0[uniform[J,L,L],uniform[K,L,L],uniform[\[CapitalGamma],L,L],h,\[Chi]v,\[Omega]v,L,L];
E\[Lambda]=EnLagMF[uniform[J,L,L],uniform[K,L,L],uniform[\[CapitalGamma],L,L],h,\[Chi]v,\[Omega]v,L,L];
];                                        (* <- EnMF0 ?  *)
Esum=Sum[  Total[Select[Eigenvalues[HMFk[J,K,\[CapitalGamma],h,\[Chi],\[Omega],\[Eta],kTable[[l]] ]  ],#<0&]] ,{l,1,Nc}]/(2Nc); 

	dataToFile[parameters[[1,p]],L,acuracy,{j,L,\[Chi],\[Omega],{0,0},{{EMF},{Esum},{EMF+E\[Lambda]},\[CapitalDelta]seq,\[CapitalDelta]\[Omega]seq}},"free",NbName]; 
Print["Max Step = ", j,"; Delta=",round\[CapitalDelta]@\[CapitalDelta]1(*,"; \[CapitalDelta]t = ",IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]*),"; E=",{EMF,Esum,EMF+E\[Lambda]},";  p=",p,"/",Length@parameters[[1]]];
{jG,LG,\[Chi]G,\[Omega]G,\[Xi]G,EnG}= loadData[toPath[parameters[[1,p]],L,acuracy,"free",NbName ]  ];
 ];  , {p,1,Length[parameters[[1]]]}  ],{l,1,Length@Ls}];Print[" "];


(* ::Subsubsection::Bold::Closed:: *)
(*vortex free + gradually increase  electric field*)


(*Module[{\[CapitalDelta]t}, t0v=AbsoluteTime[];
\[CapitalDelta]t= UnitConvert[ Quantity[N[t0v-tvf], "Seconds" ], "Minutes" ];   
Print["Free loop timing= ", IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]    ];t0=t0v;Print[" "] Print[" "];];
Print["    Starting vortex free + electric field loop: "];Print[" "]*)


(* Begin comment 
minSteps=10;
Do[   Module[{\[Chi]G,\[Omega]G,jG,LG,EnG ,gauge="g0"},   (* <-  the 1st difference : g0 \[UndirectedEdge] g4 *)
Module[{ J,K,\[CapitalGamma],Jmod,Kmod,\[CapitalGamma]mod,Jv,Kv,\[CapitalGamma]v,L=Ls[[l]],Nc,h ,\[CapitalLambda],T,En,EMF,Esum,E\[Lambda],EnList={{},{},{}},\[Xi]G,\[CapitalDelta]seq={},\[CapitalDelta]\[Omega]seq={},\[CapitalDelta]\[Omega],u2=0,u1,u0,\[Chi]={0,0,0},\[Omega]={0,0},\[Xi]={0,0},j,\[CapitalDelta]1=1,\[CapitalDelta]2=2.56,ES,gap,\[CapitalDelta]t,hp=Mod[p,Length@hV,1] }, 
{J,K,\[CapitalGamma],h,Jmod,Kmod,\[CapitalGamma]mod}=parameters[[1,p]][[1;;7]]; 
Nc=L^2;
\[Omega]=Table[\[Omega]GA,2,{r,1,Nc} ];
T=Tmat[L,L];
Kv=uniform[K,L,L];

u0=uniformU[-1,L]; (* <-  the 2nd difference : gauge4v *)

{jG,LG,\[Chi]G,\[Omega]G,\[Xi]G,EnG}= loadData[toPath[parameters[[1,p]],L,acuracy,"free",NbName ]  ]; 
(*\[Omega]G=sym\[Omega][\[Omega]G];*)
\[Omega][[1]]=Table[\[Omega]G[[1]],{r,1,Nc} ];
\[Omega][[2]]=Table[\[Omega]G[[2]],{r,1,Nc} ];
\[Chi][[1]]=Table[\[Chi]G[[1]],{r,1,Nc} ];
\[Chi][[2]]=Table[\[Chi]G[[2]],{r,1,Nc} ];
\[Chi][[3]]=Table[\[Chi]G[[3]],{r,1,Nc} ];

(* Print[" for Pure Kitaev model: " ];*)
Module[ {h0=Norm[h],\[Kappa]0,\[Kappa],\[Lambda]=to\[Lambda][h],\[Chi]0,Hpure,Tpure,Upure,Epure},   
\[Kappa]0=(*(h0/Sqrt[3])^3/(0.262)^2*) toKappa[h]; 
\[Kappa]=N@(Round[10000 \[Kappa]0]/10000);   
\[Chi]0={0,0,0};  (*
\[Omega]0=Table[\[Omega]GA,2,{r,1,Nc} ];*)
Hpure=Hreal[Kv,\[Kappa],\[Lambda],u0,L,0, {0,0}] ;   
Tpure=TmatPure[L];  
Upure=UmatPure[Tpure\[ConjugateTranspose] . Hpure . Tpure];
Epure=Total[Select[Quiet@Eigenvalues[Hpure],#<0&]]/(Nc);
{\[Chi]0[[1]],\[Chi]0[[2]],\[Chi]0[[3]]}=toMFparametersPure[Upure,u0,L,0];   
dataToFilePure[parameters[[1,p]],L,acuracy,gauge,{0,L,\[Chi]0,{{},{}},{{{},{},{}},{{},{},{}}},{Epure} } ]; 
Print["Pure data saved: "];
Print["Kappa=",\[Kappa],"; Lambda=", \[Lambda], "; Epure=",Epure,"; "];Print[];
Do[ \[Chi][[1,r]][[1,1]]=\[Chi]0[[1,r]][[1,1]]; \[Chi][[2,r]][[1,1]]=\[Chi]0[[2,r]][[1,1]]; \[Chi][[3,r]][[1,1]]=\[Chi]0[[3,r]][[1,1]]; ,{r,1,Nc}];
  ]; 

(*\[Chi]=\[Chi]gauge4v[\[Chi],L];*)   (* <-  the 3rd difference :  \[Chi]gauge4v *)

 Do[ {J,K,\[CapitalGamma],h,Jmod,Kmod,\[CapitalGamma]mod}=parameters[[ev,p]][[1;;7]]; Print[" "];
Print["J=",J, "; K=",K, "; G=",\[CapitalGamma],"; Jmod=",Jmod, "; Kmod=",Kmod, "; Gmod=",\[CapitalGamma]mod, "; L=",L, "; h=(", hV[[ hp,1 ]],",",hV[[ hp,2 ]],",",hV[[ hp,3]],"); eV=",N[Round[1000 eVs[[ev]]/1700]/1000]"; "];Print[" "];
Kv=uniform[K,L,L];Kv=add4VorticesMaxSpaced[ Kv,Kmod,L];
Jv=uniform[J,L,L];Jv=add4VorticesMaxSpaced[ Jv,Jmod,L];
\[CapitalGamma]v=uniform[\[CapitalGamma],L,L];\[CapitalGamma]v=add4VorticesMaxSpaced[ \[CapitalGamma]v,\[CapitalGamma]mod,L];

For[j=1, ( j<steps)\[And]((j<minSteps)\[Or](Chop[ \[CapitalDelta]1, 10^-acuracy ]!= 0)    ) , j++,   
If[j>=2,\[CapitalDelta]2=\[CapitalDelta]1; \[CapitalDelta]1=Max[ Abs@(u1-u2) ]; \[CapitalDelta]seq={\[CapitalDelta]seq,{j,\[CapitalDelta]1}};  \[CapitalDelta]\[Omega]=1/(2 Nc) Sum[Abs[\[Omega][[\[Sigma],r,2,1]]+\[Omega][[\[Sigma],r,3,4]]],{r,1,Nc},{\[Sigma],1,2}]; \[CapitalDelta]\[Omega]seq={\[CapitalDelta]\[Omega]seq,{j,\[CapitalDelta]\[Omega]}}    ]; 

Module[{H,u,TUh,Heff,\[Lambda]1,\[Lambda]2,loaddata},
If[j<=1, loaddata=loadDataTry[toPath[parameters[[ev,p]],L,acuracy,gauge,NbName]  ];
If[!(loaddata===$Failed),{j,L,\[Chi],\[Omega],\[Xi],{EnList[[1]],EnList[[2]],EnList[[3]],\[CapitalDelta]seq,\[CapitalDelta]\[Omega]seq}}=loaddata]];

Heff=HeffList[Jv,Kv,\[CapitalGamma]v,h,\[Omega]];(*  \[Lambda]1=1/2 Heff;\[Lambda]2=1/2 Heff;  *)
\[Lambda]1=\[Lambda]1List[Heff,\[Omega]];\[Lambda]2=\[Lambda]2List[Heff,\[Omega]]; 
H=HMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L,\[Lambda]1,\[Lambda]2,Heff]; 
u=Umat[T\[ConjugateTranspose] . H . T];
u1=Re@Chop@icc[u,L,T];

	EMF=EnMF0[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];                                        (* <- EnMF0 ?  *)
	Esum=Total[Select[Quiet@Eigenvalues[H],#<0&]]/(2Nc);
	E\[Lambda]=EnLagMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];
	EnList[[1]]={EnList[[1]],{j,EMF}};
    EnList[[2]]={EnList[[2]],{j,Esum}};
    EnList[[3]]={EnList[[3]],{j,EMF+E\[Lambda]}};   
Do[ Module[{m,n,\[Alpha],\[Beta],rz,rx,ry,Io}, rz=\[LeftFloor]R0/16\[RightFloor];\[Beta]=\[LeftFloor]1+(R0-16rz)/4\[RightFloor];\[Alpha]=R0-16rz-4(\[Beta]-1)+1;
   n=\[LeftFloor]rz/L\[RightFloor];m=rz-n L;rx=Mod[m+1,L]+n L;ry=m+Mod[n+1,L] L;Io=Mod[\[Alpha]+8rz,8Nc,1]; 
    \[Chi][[1,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rx,8Nc,1],Io]];
    \[Chi][[2,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8ry,8Nc,1],Io]];
    \[Chi][[3,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rz,8Nc,1],Io]];
    \[Omega][[1,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+8rz,8Nc,1],Io]];
    \[Omega][[2,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rz,8Nc,1],Mod[Io+4,8Nc,1] ]];
    ]; , {R0,0,16Nc-1}  ];    
]; 
u2=u1;        
Print[" j =",j, "/",steps, "; Delta=",round\[CapitalDelta]@\[CapitalDelta]1";  E=", N[Round[1000000(#)]/1000000]&@{EMF,Esum,EMF+E\[Lambda]},"; "  ];     
         ];  
Module[{H,u,Heff,\[Lambda]1,\[Lambda]2},Heff=HeffList[Jv,Kv,\[CapitalGamma]v,h,\[Omega]]; (*
\[Lambda]1=1/2 Heff;\[Lambda]2=1/2 Heff;*)
\[Lambda]1=\[Lambda]1List[Heff,\[Omega]];\[Lambda]2=\[Lambda]2List[Heff,\[Omega]]; 
H=HMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L,\[Lambda]1,\[Lambda]2,Heff]; 
    u=Umat[T\[ConjugateTranspose] . H . T];
	{\[Chi][[1]],\[Chi][[2]],\[Chi][[3]],\[Omega][[1]],\[Omega][[2]],\[Xi][[1]],\[Xi][[2]],u1}=BiParallel[u,L,T];

	EMF=EnMF0[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];                                        (* <- EnMF0 ?  *)
	Esum=Total[Select[Quiet@Eigenvalues[H],#<0&]]/(2Nc);
	E\[Lambda]=EnLagMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];
	EnList[[1]]={EnList[[1]],{j,EMF}};
    EnList[[2]]={EnList[[2]],{j,Esum}};
    EnList[[3]]={EnList[[3]],{j,EMF+E\[Lambda]}};      

];  	\[CapitalDelta]seq=Partition[Flatten[\[CapitalDelta]seq],2];
	EnList[[1]]=Partition[Flatten[EnList[[1]]],2];
	EnList[[2]]=Partition[Flatten[EnList[[2]]],2];
	EnList[[3]]=Partition[Flatten[EnList[[3]]],2]; 
	\[CapitalDelta]\[Omega]seq=Partition[Flatten[\[CapitalDelta]\[Omega]seq],2]; 
	dataToFile[parameters[[ev,p]],L,acuracy,{j,L,\[Chi],\[Omega],\[Xi],{EnList[[1]],EnList[[2]],EnList[[3]],\[CapitalDelta]seq,\[CapitalDelta]\[Omega]seq}},gauge,NbName];     
Print[ "ev=",ev ,"/", Length@eVs"; j MAX=",j, "/",steps, "; Delta=",\[CapitalDelta]1, "; "   ];
, {ev,1,Length[parameters]} ]   ;


	t1=AbsoluteTime[];\[CapitalDelta]t= UnitConvert[ Quantity[N[t1 -t0], "Seconds" ], "Hours" ];t0=t1;Print[" "];
	Print[ "p=",p,"/",Length@parameters[[1]], "; l=",l, "/",Length@Ls, "; \[CapitalDelta]t = ",IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ]   ]


  ];
]  , {l,1,Length@Ls}, {p,1,Length[parameters[[1]] ]}  ];      
                                        
end comment *) 


(* ::Subsubsection::Bold::Closed:: *)
(*four vortex + gradually increase  electric field*)


(*Module[{\[CapitalDelta]t},t4v=AbsoluteTime[];
\[CapitalDelta]t= UnitConvert[ Quantity[N[t4v-t0v], "Seconds" ], "Minutes" ];   
Print["Free loop + electric field timing= ", IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]    ];t0=t4v;Print[" "] Print[" "];];
Print["    Starting four vortex + electric field loop: "];Print[" "]*)


(* ::Code::Bold:: *)
(**)


(*minSteps=10;
Do[   Module[{\[Chi]G,\[Omega]G,jG,LG,EnG ,gauge="g4"},   (* <-  the 1st difference : g0 \[UndirectedEdge] g4 *)
Module[{ J,K,\[CapitalGamma],Jmod,Kmod,\[CapitalGamma]mod,Jv,Kv,\[CapitalGamma]v,L=Ls[[l]],Nc,h ,\[CapitalLambda],T,En,EMF,E\[Lambda],EnList={{},{},{}},\[Xi]G,\[CapitalDelta]seq={},\[CapitalDelta]\[Omega]seq={},\[CapitalDelta]\[Omega],u2,u1,u0,\[Chi]={0,0,0},\[Omega]={0,0},\[Xi]={0,0},j,\[CapitalDelta]1=1,\[CapitalDelta]2=2.56,ES,gap,\[CapitalDelta]t,hp=Mod[p,Length@hV,1] }, 
{J,K,\[CapitalGamma],h,Jmod,Kmod,\[CapitalGamma]mod}=parameters[[1,p]][[1;;7]];  
Nc=L^2;
\[Omega]=Table[\[Omega]GA,2,{r,1,Nc} ];
T=Tmat[L,L];
Kv=uniform[K,L,L];

u0=gauge4v[uniformU[-1,L],L]; (* <-  the 2nd difference : gauge4v *)

{jG,LG,\[Chi]G,\[Omega]G,\[Xi]G,EnG}= loadData[toPath[parameters[[1,p]],L,acuracy,"free",NbName ]  ]; 
(*\[Omega]G=sym\[Omega][\[Omega]G];*)
\[Omega][[1]]=Table[\[Omega]G[[1]],{r,1,Nc} ];
\[Omega][[2]]=Table[\[Omega]G[[2]],{r,1,Nc} ];
\[Chi][[1]]=Table[\[Chi]G[[1]],{r,1,Nc} ];
\[Chi][[2]]=Table[\[Chi]G[[2]],{r,1,Nc} ];
\[Chi][[3]]=Table[\[Chi]G[[3]],{r,1,Nc} ];

(* Print[" for Pure Kitaev model: " ];*)
Module[ {h0=Norm[h],\[Kappa]0,\[Kappa],\[Chi]0,Hpure,Tpure,Upure},   
\[Kappa]0=(*(h0/Sqrt[3])^3/(0.262)^2*) toKappa[h]; 
\[Kappa]=N@(Round[10000 \[Kappa]0]/10000);   
\[Chi]0={0,0,0};  (*
\[Omega]0=Table[\[Omega]GA,2,{r,1,Nc} ];*)
Hpure=Hreal[Kv,\[Kappa],to\[Lambda][h],u0,L,0, {0,0}] ;   
Tpure=TmatPure[L];  
Upure=UmatPure[Tpure\[ConjugateTranspose] . Hpure . Tpure];
Epure=Total[Select[Quiet@Eigenvalues[Hpure],#<0&]]/(Nc);
{\[Chi]0[[1]],\[Chi]0[[2]],\[Chi]0[[3]]}=toMFparametersPure[Upure,u0,L,0];   
dataToFilePure[parameters[[ev,p]],L,acuracy,gauge,{0,L,\[Chi]0,{{},{}},{{{},{},{}},{{},{},{}}},{Epure} } ]; 
Print["Pure data saved: "];
Print["Kappa=",\[Kappa]," Lambda=", \[Lambda]"; Epure=",Epure,"; "];Print[];
Do[ \[Chi][[1,r]][[1,1]]=\[Chi]0[[1,r]][[1,1]]; \[Chi][[2,r]][[1,1]]=\[Chi]0[[2,r]][[1,1]]; \[Chi][[3,r]][[1,1]]=\[Chi]0[[3,r]][[1,1]]; ,{r,1,Nc}];
  ]; 

\[Chi]=\[Chi]gauge4v[\[Chi],L];   (* <-  the 3rd difference :  \[Chi]gauge4v *)

 Do[ {J,K,\[CapitalGamma],h,Jmod,Kmod,\[CapitalGamma]mod}=parameters[[ev,p]][[1;;7]]; Print[" "];
Print["J=",J, "; K=",K, "; G=",\[CapitalGamma],"; Jmod=",Jmod, "; Kmod=",Kmod, "; Gmod=",\[CapitalGamma]mod, "; L=",L, "; h=(", hV[[ hp,1 ]],",",hV[[ hp,2 ]],",",hV[[ hp,3]],"); eV=",N[Round[1000 eVs[[ev]]/1700]/1000]"; "];Print[" "];
Kv=uniform[K,L,L];Kv=add4VorticesMaxSpaced[ Kv,Kmod,L];
Jv=uniform[J,L,L];Jv=add4VorticesMaxSpaced[ Jv,Jmod,L];
\[CapitalGamma]v=uniform[\[CapitalGamma],L,L];\[CapitalGamma]v=add4VorticesMaxSpaced[ \[CapitalGamma]v,\[CapitalGamma]mod,L];
For[j=1, ( j<steps)\[And]((j<minSteps)\[Or](Chop[ \[CapitalDelta]1, 10^-acuracy ]!= 0)    ) , j++,   

Module[{H,u,TUh,Heff,\[Lambda]1,\[Lambda]2},Heff=HeffList[Jv,Kv,\[CapitalGamma]v,h,\[Omega]];
\[Lambda]1=1/2 Heff;\[Lambda]2=1/2 Heff;
(*\[Lambda]1=\[Lambda]1List[Heff,\[Omega]];\[Lambda]2=\[Lambda]2List[Heff,\[Omega]]; *)
H=HMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L,\[Lambda]1,\[Lambda]2,Heff]; 
u=Umat[T\[ConjugateTranspose] . H . T];
u1=Re@Chop@icc[u,L,T];

	EMF=EnMF0[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];                                        (* <- EnMF0 ?  *)
	E\[Lambda]=EnLagMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];
	EnList[[1]]={EnList[[1]],{j,EMF}};
    EnList[[2]]={EnList[[2]],{j,0}};
    EnList[[3]]={EnList[[3]],{j,EMF+E\[Lambda]}};   

Do[ Module[{m,n,\[Alpha],\[Beta],rz,rx,ry,Io}, rz=\[LeftFloor]R0/16\[RightFloor];\[Beta]=\[LeftFloor]1+(R0-16rz)/4\[RightFloor];\[Alpha]=R0-16rz-4(\[Beta]-1)+1;
   n=\[LeftFloor]rz/L\[RightFloor];m=rz-n L;rx=Mod[m+1,L]+n L;ry=m+Mod[n+1,L] L;Io=Mod[\[Alpha]+8rz,8Nc,1]; 
    \[Chi][[1,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rx,8Nc,1],Io]];
    \[Chi][[2,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8ry,8Nc,1],Io]];
    \[Chi][[3,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rz,8Nc,1],Io]];
    \[Omega][[1,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+8rz,8Nc,1],Io]];
    \[Omega][[2,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rz,8Nc,1],Mod[Io+4,8Nc,1] ]];
    ]; , {R0,0,16Nc-1}  ];     
]; 
If[j>=2,\[CapitalDelta]2=\[CapitalDelta]1; \[CapitalDelta]1=Max[ Abs@(u1-u2) ]; \[CapitalDelta]seq={\[CapitalDelta]seq,{j,\[CapitalDelta]1}};  \[CapitalDelta]\[Omega]=1/(2 Nc) Sum[Abs[\[Omega][[\[Sigma],r,2,1]]+\[Omega][[\[Sigma],r,3,4]]],{r,1,Nc},{\[Sigma],1,2}]; \[CapitalDelta]\[Omega]seq={\[CapitalDelta]\[Omega]seq,{j,\[CapitalDelta]\[Omega]}}    ]; 
u2=u1;        
Print[" j =",j, "/",steps, "; Delta=",round\[CapitalDelta]@\[CapitalDelta]1";  E=", N[Round[1000000(EMF)]/1000000],"; "  ];     
         ];  
Module[{H,u,Heff,\[Lambda]1,\[Lambda]2},Heff=HeffList[Jv,Kv,\[CapitalGamma]v,h,\[Omega]]; (*
\[Lambda]1=1/2 Heff;\[Lambda]2=1/2 Heff;*)
\[Lambda]1=\[Lambda]1List[Heff,\[Omega]];\[Lambda]2=\[Lambda]2List[Heff,\[Omega]]; 
H=HMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L,\[Lambda]1,\[Lambda]2,Heff]; 
    u=Umat[T\[ConjugateTranspose] . H . T];
	{\[Chi][[1]],\[Chi][[2]],\[Chi][[3]],\[Omega][[1]],\[Omega][[2]],\[Xi][[1]],\[Xi][[2]],u1}=BiParallel[u,L,T]; 
	EMF=EnMF0[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L]; 
	E\[Lambda]=EnLagMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];
	EnList[[1]]={EnList[[1]],{j,EMF}};
    EnList[[2]]={EnList[[2]],{j,0}};
    EnList[[3]]={EnList[[3]],{j,EMF+E\[Lambda]}};           
];  	\[CapitalDelta]seq=Partition[Flatten[\[CapitalDelta]seq],2];
	EnList[[1]]=Partition[Flatten[EnList[[1]]],2];
	EnList[[2]]=Partition[Flatten[EnList[[2]]],2];
	EnList[[3]]=Partition[Flatten[EnList[[3]]],2]; 
	\[CapitalDelta]\[Omega]seq=Partition[Flatten[\[CapitalDelta]\[Omega]seq],2]; 
	dataToFile[parameters[[ev,p]],L,acuracy,{j,L,\[Chi],\[Omega],\[Xi],{EnList[[1]],EnList[[2]],EnList[[3]],\[CapitalDelta]seq,\[CapitalDelta]\[Omega]seq}},gauge,NbName];     
Print[ "ev=",ev ,"/", Length@eVs"; j MAX=",j, "/",steps, "; Delta=",\[CapitalDelta]1, "; "   ]
, {ev,1,Length[parameters]} ]   ;


	t1=AbsoluteTime[];\[CapitalDelta]t= UnitConvert[ Quantity[N[t1 -t0], "Seconds" ], "Hours" ];t0=t1; Print[" "];
	Print[ "p=",p,"/",Length@parameters[[1]], "; l=",l, "/",Length@Ls, "; \[CapitalDelta]t = ",IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ]   ]



  ];
]  , {l,1,Length@Ls}, {p,1,Length[parameters[[1]] ]}  ]                                                 *)


(*Module[{\[CapitalDelta]t},t1=AbsoluteTime[];\[CapitalDelta]t= UnitConvert[ Quantity[N[t1-t4v], "Seconds" ], "Hours" ];
Print[ "4 vortices loop timing \[CapitalDelta]t = ",IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ]   ] ];*)


(* ::Subsubsection::Bold:: *)
(*vortex free + gradually increase  parameters *)


Module[{\[CapitalDelta]t}, t0v=AbsoluteTime[];
\[CapitalDelta]t= UnitConvert[ Quantity[N[t0v-tvf], "Seconds" ], "Minutes" ];   
Print["Free loop timing= ", IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]    ];t0=t0v;Print[" "] Print[" "];];
Print["    Starting vortex free + electric field loop: "];Print[" "]


(*minSteps=10;
Do[   Module[{\[Chi]G,\[Omega]G,jG,LG,EnG ,gauge="g0"},   (* <-  the 1st difference : g0 \[UndirectedEdge] g4 *)
Module[{ J,K,\[CapitalGamma],Jmod,Kmod,\[CapitalGamma]mod,Jv,Kv,\[CapitalGamma]v,L=Ls[[l]],Nc,h ,\[CapitalLambda],T,En,EMF,Esum,E\[Lambda],EnList={{},{},{}},\[Xi]G,\[CapitalDelta]seq={},\[CapitalDelta]\[Omega]seq={},\[CapitalDelta]\[Omega],u2,u1,u0,\[Chi]={0,0,0},\[Omega]={0,0},\[Xi]={0,0},j,\[CapitalDelta]1=1,\[CapitalDelta]2=2.56,ES,gap,\[CapitalDelta]t,hp}, 
{J,K,\[CapitalGamma],h,Jmod,Kmod,\[CapitalGamma]mod}=parameters[[ev,1]][[1;;7]];  
Nc=L^2;
\[Omega]=Table[\[Omega]GA,2,{r,1,Nc} ];
T=Tmat[L,L];
Kv=uniform[K,L,L];

u0=uniformU[-1,L]; (* <-  the 2nd difference : gauge4v *)

{jG,LG,\[Chi]G,\[Omega]G,\[Xi]G,EnG}= loadData[toPath[parameters[[1,1]],L,acuracy,"free",NbName ]  ]; 
(*\[Omega]G=sym\[Omega][\[Omega]G];*)
\[Omega][[1]]=Table[\[Omega]G[[1]],{r,1,Nc} ];
\[Omega][[2]]=Table[\[Omega]G[[2]],{r,1,Nc} ];
\[Chi][[1]]=Table[\[Chi]G[[1]],{r,1,Nc} ];
\[Chi][[2]]=Table[\[Chi]G[[2]],{r,1,Nc} ];
\[Chi][[3]]=Table[\[Chi]G[[3]],{r,1,Nc} ];

(* Print[" for Pure Kitaev model: " ];*)
Module[ {h0=Norm[h],\[Kappa]0,\[Kappa],\[Lambda]=to\[Lambda][h],\[Chi]0,Hpure,Tpure,Upure,Epure},   
\[Kappa]0=(*(h0/Sqrt[3])^3/(0.262)^2*) toKappa[h]; 
\[Kappa]=N@(Round[10000 \[Kappa]0]/10000);   
\[Chi]0={0,0,0};  (*
\[Omega]0=Table[\[Omega]GA,2,{r,1,Nc} ];*)
Hpure=Hreal[Kv,\[Kappa],\[Lambda],u0,L,0, {0,0}] ;   
Tpure=TmatPure[L];  
Upure=UmatPure[Tpure\[ConjugateTranspose] . Hpure . Tpure];
Epure=Total[Select[Quiet@Eigenvalues[Hpure],#<0&]]/(Nc);
{\[Chi]0[[1]],\[Chi]0[[2]],\[Chi]0[[3]]}=toMFparametersPure[Upure,u0,L,0];   
dataToFilePure[parameters[[ev,1]],L,acuracy,gauge,{0,L,\[Chi]0,{{},{}},{{{},{},{}},{{},{},{}}},{Epure} } ]; 
Print["Pure data saved: "];
Print["Kappa=",\[Kappa],"; Lambda=", \[Lambda] ,"; Epure=",Epure,"; "];Print[];
Do[ \[Chi][[1,r]][[1,1]]=\[Chi]0[[1,r]][[1,1]]; \[Chi][[2,r]][[1,1]]=\[Chi]0[[2,r]][[1,1]]; \[Chi][[3,r]][[1,1]]=\[Chi]0[[3,r]][[1,1]]; ,{r,1,Nc}];
  ]; 
 (* <-  the 3rd difference :  \[Chi]gauge4v *)

 Do[ {J,K,\[CapitalGamma],h,Jmod,Kmod,\[CapitalGamma]mod}=parameters[[ev,p]][[1;;7]]; Print[" "];   hp=Mod[p,Length@hV,1] ;
Print["J=",J, "; K=",K, "; G=",\[CapitalGamma],"; Jmod=",Jmod, "; Kmod=",Kmod, "; Gmod=",\[CapitalGamma]mod, "; L=",L, "; h=(", hV[[ hp,1 ]],",",hV[[ hp,2 ]],",",hV[[ hp,3]],"); eV=",N[Round[1000 eVs[[ev]]/1700]/1000]"; "];Print[" "];
Kv=uniform[K,L,L];Kv=add4VorticesMaxSpaced[ Kv,Kmod,L];
Jv=uniform[J,L,L];Jv=add4VorticesMaxSpaced[ Jv,Jmod,L];
\[CapitalGamma]v=uniform[\[CapitalGamma],L,L];\[CapitalGamma]v=add4VorticesMaxSpaced[ \[CapitalGamma]v,\[CapitalGamma]mod,L];

For[j=1, ( j<steps)\[And]((j<minSteps)\[Or](Chop[ \[CapitalDelta]1, 10^-acuracy ]!= 0)    ) , j++,   

Module[{H,u,TUh,Heff,\[Lambda]1,\[Lambda]2,loaddata},
(*If[j<=1, loaddata=loadDataTry[toPath[parameters[[ev,p]],L,acuracy,gauge,NbName]  ];
If[!(loaddata===$Failed),{j,L,\[Chi],\[Omega],\[Xi],{EnList[[1]],EnList[[2]],EnList[[3]],\[CapitalDelta]seq,\[CapitalDelta]\[Omega]seq}}=loaddata; Continue[]] ];*)

Heff=HeffList[Jv,Kv,\[CapitalGamma]v,h,\[Omega]];
\[Lambda]1=1/2 Heff;\[Lambda]2=1/2 Heff;
(*\[Lambda]1=\[Lambda]1List[Heff,\[Omega]];\[Lambda]2=\[Lambda]2List[Heff,\[Omega]]; *)
H=HMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L,\[Lambda]1,\[Lambda]2,Heff]; 
u=Umat[T\[ConjugateTranspose] . H . T];
u1=Re@Chop@icc[u,L,T];

	EMF=EnMF0[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];                                        (* <- EnMF0 ?  *)
	Esum=Total[Select[Quiet@Eigenvalues[H],#<0&]]/(2Nc);
	E\[Lambda]=EnLagMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];
	EnList[[1]]={EnList[[1]],{j,EMF}};
    EnList[[2]]={EnList[[2]],{j,Esum}};
    EnList[[3]]={EnList[[3]],{j,EMF+E\[Lambda]}};   

Do[ Module[{m,n,\[Alpha],\[Beta],rz,rx,ry,Io}, rz=\[LeftFloor]R0/16\[RightFloor];\[Beta]=\[LeftFloor]1+(R0-16rz)/4\[RightFloor];\[Alpha]=R0-16rz-4(\[Beta]-1)+1;
   n=\[LeftFloor]rz/L\[RightFloor];m=rz-n L;rx=Mod[m+1,L]+n L;ry=m+Mod[n+1,L] L;Io=Mod[\[Alpha]+8rz,8Nc,1]; 
    \[Chi][[1,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rx,8Nc,1],Io]];
    \[Chi][[2,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8ry,8Nc,1],Io]];
    \[Chi][[3,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rz,8Nc,1],Io]];
    \[Omega][[1,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+8rz,8Nc,1],Io]];
    \[Omega][[2,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rz,8Nc,1],Mod[Io+4,8Nc,1] ]];
    ]; , {R0,0,16Nc-1}  ];     
]; 
If[j>=2,\[CapitalDelta]2=\[CapitalDelta]1; \[CapitalDelta]1=Max[ Abs@(u1-u2) ]; \[CapitalDelta]seq={\[CapitalDelta]seq,{j,\[CapitalDelta]1}};  \[CapitalDelta]\[Omega]=1/(2 Nc) Sum[Abs[\[Omega][[\[Sigma],r,2,1]]+\[Omega][[\[Sigma],r,3,4]]],{r,1,Nc},{\[Sigma],1,2}]; \[CapitalDelta]\[Omega]seq={\[CapitalDelta]\[Omega]seq,{j,\[CapitalDelta]\[Omega]}}    ]; 
u2=u1;        
Print[" j =",j, "/",steps, "; Delta=",round\[CapitalDelta]@\[CapitalDelta]1";  E=", N[Round[1000000(#)]/1000000]&@{EMF,Esum,EMF+E\[Lambda]},"; "  ];     
         ];  
Module[{H,u,Heff,\[Lambda]1,\[Lambda]2},Heff=HeffList[Jv,Kv,\[CapitalGamma]v,h,\[Omega]]; (*
\[Lambda]1=1/2 Heff;\[Lambda]2=1/2 Heff;*)
\[Lambda]1=\[Lambda]1List[Heff,\[Omega]];\[Lambda]2=\[Lambda]2List[Heff,\[Omega]]; 
H=HMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L,\[Lambda]1,\[Lambda]2,Heff]; 
    u=Umat[T\[ConjugateTranspose] . H . T];
	{\[Chi][[1]],\[Chi][[2]],\[Chi][[3]],\[Omega][[1]],\[Omega][[2]],\[Xi][[1]],\[Xi][[2]],u1}=BiParallel[u,L,T]; 
	EMF=EnMF0[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L]; 
	E\[Lambda]=EnLagMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];
	Esum=Total[Select[Quiet@Eigenvalues[H],#<0&]]/(2Nc);
	EnList[[1]]={EnList[[1]],{j,EMF}};
    EnList[[2]]={EnList[[2]],{j,Esum}};
    EnList[[3]]={EnList[[3]],{j,EMF+E\[Lambda]}};           
];  	\[CapitalDelta]seq=Partition[Flatten[\[CapitalDelta]seq],2];
	EnList[[1]]=Partition[Flatten[EnList[[1]]],2];
	EnList[[2]]=Partition[Flatten[EnList[[2]]],2];
	EnList[[3]]=Partition[Flatten[EnList[[3]]],2]; 
	\[CapitalDelta]\[Omega]seq=Partition[Flatten[\[CapitalDelta]\[Omega]seq],2]; 
	dataToFile[parameters[[ev,p]],L,acuracy,{j,L,\[Chi],\[Omega],\[Xi],{EnList[[1]],EnList[[2]],EnList[[3]],\[CapitalDelta]seq,\[CapitalDelta]\[Omega]seq}},gauge,NbName];   
  
Print[ "ev=",ev ,"/", Length@eVs"; j MAX=",j, "/",steps, "; Delta=",\[CapitalDelta]1, "; "   ];
t1=AbsoluteTime[];\[CapitalDelta]t= UnitConvert[ Quantity[N[t1 -t0], "Seconds" ], "Hours" ];t0=t1; Print[" "];
Print[ "p=",p,"/",Length@parameters[[1]], "; l=",l, "/",Length@Ls, "; \[CapitalDelta]t = ",IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ]   ]


 , {p,1,Length[parameters[[1]] ]}   ];


  ];
]  , {l,1,Length@Ls} , {ev,1,Length[parameters]} ];  *)                                                 


(* ::Subsubsection::Bold:: *)
(*four vortex + gradually increase  parameters *)


Module[{\[CapitalDelta]t},t4v=AbsoluteTime[];
\[CapitalDelta]t= UnitConvert[ Quantity[N[t4v-t0v], "Seconds" ], "Minutes" ];   
Print["Free loop + electric field timing= ", IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Seconds" ]    ];t0=t4v;Print[" "] Print[" "];];
Print["    Starting four vortex -fixed eV- gradually changing parameters  "];Print[" "]


minSteps=10;
Do[   Module[{\[Chi]G,\[Omega]G,jG,LG,EnG ,gauge="g4"},   (* <-  the 1st difference : g0 \[UndirectedEdge] g4 *)
Module[{ J,K,\[CapitalGamma],Jmod,Kmod,\[CapitalGamma]mod,Jv,Kv,\[CapitalGamma]v,L=Ls[[l]],Nc,h ,\[CapitalLambda],T,En,EMF,Esum,E\[Lambda],EnList={{},{},{}},\[Xi]G,\[CapitalDelta]seq={},\[CapitalDelta]\[Omega]seq={},\[CapitalDelta]\[Omega],u2,u1,u0,\[Chi]={0,0,0},\[Omega]={0,0},\[Xi]={0,0},j,\[CapitalDelta]1=1,\[CapitalDelta]2=2.56,ES,gap,\[CapitalDelta]t,hp}, 
{J,K,\[CapitalGamma],h,Jmod,Kmod,\[CapitalGamma]mod}=parameters[[ev,1]][[1;;7]];  
Nc=L^2;
\[Omega]=Table[\[Omega]GA,2,{r,1,Nc} ];
T=Tmat[L,L];
Kv=uniform[K,L,L];

u0=gauge4v[uniformU[-1,L],L]; (* <-  the 2nd difference : gauge4v *)

{jG,LG,\[Chi]G,\[Omega]G,\[Xi]G,EnG}= loadData[toPath[parameters[[1,1]],L,acuracy,"free",NbName ]  ]; 
(*\[Omega]G=sym\[Omega][\[Omega]G];*)
\[Omega][[1]]=Table[\[Omega]G[[1]],{r,1,Nc} ];
\[Omega][[2]]=Table[\[Omega]G[[2]],{r,1,Nc} ];
\[Chi][[1]]=Table[\[Chi]G[[1]],{r,1,Nc} ];
\[Chi][[2]]=Table[\[Chi]G[[2]],{r,1,Nc} ];
\[Chi][[3]]=Table[\[Chi]G[[3]],{r,1,Nc} ];

(* Print[" for Pure Kitaev model: " ];*)
Module[ {h0=Norm[h],\[Kappa]0,\[Kappa],\[Lambda]=to\[Lambda][h],\[Chi]0,Hpure,Tpure,Upure,Epure},   
\[Kappa]0=(*(h0/Sqrt[3])^3/(0.262)^2*) toKappa[h]; 
\[Kappa]=N@(Round[10000 \[Kappa]0]/10000);   
\[Chi]0={0,0,0};  (*
\[Omega]0=Table[\[Omega]GA,2,{r,1,Nc} ];*)
Hpure=Hreal[Kv,\[Kappa],\[Lambda],u0,L,0, {0,0}] ;   
Tpure=TmatPure[L];  
Upure=UmatPure[Tpure\[ConjugateTranspose] . Hpure . Tpure];
Epure=Total[Select[Quiet@Eigenvalues[Hpure],#<0&]]/(Nc);
{\[Chi]0[[1]],\[Chi]0[[2]],\[Chi]0[[3]]}=toMFparametersPure[Upure,u0,L,0];   
dataToFilePure[parameters[[ev,1]],L,acuracy,gauge,{0,L,\[Chi]0,{{},{}},{{{},{},{}},{{},{},{}}},{Epure} } ]; 
Print["Pure data saved: "];
Print["Kappa=",\[Kappa],"; Lambda=", \[Lambda] ,"; Epure=",Epure,"; "];Print[];
Do[ \[Chi][[1,r]][[1,1]]=\[Chi]0[[1,r]][[1,1]]; \[Chi][[2,r]][[1,1]]=\[Chi]0[[2,r]][[1,1]]; \[Chi][[3,r]][[1,1]]=\[Chi]0[[3,r]][[1,1]]; ,{r,1,Nc}];
  ]; 
\[Chi]=\[Chi]gauge4v[\[Chi],L];   (* <-  the 3rd difference :  \[Chi]gauge4v *)

 Do[  p+=1;

{J,K,\[CapitalGamma],h,Jmod,Kmod,\[CapitalGamma]mod}=parameters[[ev,p]][[1;;7]]; Print[" "];   hp=Mod[p,Length@hV,1] ;
Print["J=",J, "; K=",K, "; G=",\[CapitalGamma],"; Jmod=",Jmod, "; Kmod=",Kmod, "; Gmod=",\[CapitalGamma]mod, "; L=",L, "; h=(", hV[[ hp,1 ]],",",hV[[ hp,2 ]],",",hV[[ hp,3]],"); eV=",N[Round[1000 eVs[[ev]]/1700]/1000]"; "];Print[" "];
Kv=uniform[K,L,L];Kv=add4VorticesMaxSpaced[ Kv,Kmod,L];
Jv=uniform[J,L,L];Jv=add4VorticesMaxSpaced[ Jv,Jmod,L];
\[CapitalGamma]v=uniform[\[CapitalGamma],L,L];\[CapitalGamma]v=add4VorticesMaxSpaced[ \[CapitalGamma]v,\[CapitalGamma]mod,L];

For[j=1, ( j<steps)\[And]((j<minSteps)\[Or](Chop[ \[CapitalDelta]1, 10^-acuracy ]!= 0)    ) , j++,   
Print["j1=",j];
Module[{H,u,TUh,Heff,\[Lambda]1,\[Lambda]2,loaddata},
If[j<=1,   
    loaddata=loadDataTry[toPath[parameters[[ev,1]],L,acuracy,gauge,NbName] ];
    If[!(loaddata===$Failed),{j0,L0,\[Chi],\[Omega],\[Xi]0,EnG0}=loaddata;];
];

Heff=HeffList[Jv,Kv,\[CapitalGamma]v,h,\[Omega]];
\[Lambda]1=1/2 Heff;\[Lambda]2=1/2 Heff;
(*\[Lambda]1=\[Lambda]1List[Heff,\[Omega]];\[Lambda]2=\[Lambda]2List[Heff,\[Omega]]; *)
H=HMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L,\[Lambda]1,\[Lambda]2,Heff]; 
u=Umat[T\[ConjugateTranspose] . H . T];
u1=Re@Chop@icc[u,L,T];

	EMF=EnMF0[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];                                        (* <- EnMF0 ?  *)
	Esum=Total[Select[Quiet@Eigenvalues[H],#<0&]]/(2Nc);
	E\[Lambda]=EnLagMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];
	EnList[[1]]={EnList[[1]],{j,EMF}};
    EnList[[2]]={EnList[[2]],{j,Esum}};
    EnList[[3]]={EnList[[3]],{j,EMF+E\[Lambda]}};   
Print["j2=",j];

Do[ Module[{m,n,\[Alpha],\[Beta],rz,rx,ry,Io}, rz=\[LeftFloor]R0/16\[RightFloor];\[Beta]=\[LeftFloor]1+(R0-16rz)/4\[RightFloor];\[Alpha]=R0-16rz-4(\[Beta]-1)+1;
   n=\[LeftFloor]rz/L\[RightFloor];m=rz-n L;rx=Mod[m+1,L]+n L;ry=m+Mod[n+1,L] L;Io=Mod[\[Alpha]+8rz,8Nc,1]; 
    \[Chi][[1,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rx,8Nc,1],Io]];
    \[Chi][[2,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8ry,8Nc,1],Io]];
    \[Chi][[3,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rz,8Nc,1],Io]];
    \[Omega][[1,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+8rz,8Nc,1],Io]];
    \[Omega][[2,rz+1,\[Alpha],\[Beta]]]=u1[[Mod[\[Beta]+4+8rz,8Nc,1],Mod[Io+4,8Nc,1] ]];
    ]; , {R0,0,16Nc-1}  ];   
Print["j3=",j];  
]; 
If[j>=2,\[CapitalDelta]2=\[CapitalDelta]1; \[CapitalDelta]1=Max[ Abs@(u1-u2) ]; \[CapitalDelta]seq={\[CapitalDelta]seq,{j,\[CapitalDelta]1}};  \[CapitalDelta]\[Omega]=1/(2 Nc) Sum[Abs[\[Omega][[\[Sigma],r,2,1]]+\[Omega][[\[Sigma],r,3,4]]],{r,1,Nc},{\[Sigma],1,2}]; \[CapitalDelta]\[Omega]seq={\[CapitalDelta]\[Omega]seq,{j,\[CapitalDelta]\[Omega]}}    ]; 
u2=u1;        
Print[" j =",j, "/",steps, "; Delta=",round\[CapitalDelta]@\[CapitalDelta]1";  E=", N[Round[1000000(#)]/1000000]&@{EMF,Esum,EMF+E\[Lambda]},"; "  ];     
         ];  
Module[{H,u,Heff,\[Lambda]1,\[Lambda]2},Heff=HeffList[Jv,Kv,\[CapitalGamma]v,h,\[Omega]]; (*
\[Lambda]1=1/2 Heff;\[Lambda]2=1/2 Heff;*)
\[Lambda]1=\[Lambda]1List[Heff,\[Omega]];\[Lambda]2=\[Lambda]2List[Heff,\[Omega]]; 
H=HMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L,\[Lambda]1,\[Lambda]2,Heff]; 
    u=Umat[T\[ConjugateTranspose] . H . T];
	{\[Chi][[1]],\[Chi][[2]],\[Chi][[3]],\[Omega][[1]],\[Omega][[2]],\[Xi][[1]],\[Xi][[2]],u1}=BiParallel[u,L,T]; 
	EMF=EnMF0[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L]; 
	E\[Lambda]=EnLagMF[Jv,Kv,\[CapitalGamma]v,h,\[Chi],\[Omega],L,L];
	Esum=Total[Select[Quiet@Eigenvalues[H],#<0&]]/(2Nc);
	EnList[[1]]={EnList[[1]],{j,EMF}};
    EnList[[2]]={EnList[[2]],{j,Esum}};
    EnList[[3]]={EnList[[3]],{j,EMF+E\[Lambda]}};           
];  	\[CapitalDelta]seq=Partition[Flatten[\[CapitalDelta]seq],2];
	EnList[[1]]=Partition[Flatten[EnList[[1]]],2];
	EnList[[2]]=Partition[Flatten[EnList[[2]]],2];
	EnList[[3]]=Partition[Flatten[EnList[[3]]],2]; 
	\[CapitalDelta]\[Omega]seq=Partition[Flatten[\[CapitalDelta]\[Omega]seq],2]; 
	dataToFile[parameters[[ev,p]],L,acuracy,{j,L,\[Chi],\[Omega],\[Xi],{EnList[[1]],EnList[[2]],EnList[[3]],\[CapitalDelta]seq,\[CapitalDelta]\[Omega]seq}},gauge,NbName];   
  
Print[ "ev=",ev ,"/", Length@eVs"; j MAX=",j, "/",steps, "; Delta=",\[CapitalDelta]1, "; "   ];
t1=AbsoluteTime[];\[CapitalDelta]t= UnitConvert[ Quantity[N[t1 -t0], "Seconds" ], "Hours" ];t0=t1; Print[" "];
Print[ "p=",p,"/",Length@parameters[[1]], "; l=",l, "/",Length@Ls, "; \[CapitalDelta]t = ",IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ]   ]

 , {p,1,Length[parameters[[1]] ]}   ];

  ];         ]  , {l,1,Length@Ls} , {ev,1,Length[parameters]} ]                                                 


CloseKernels[];


Module[{\[CapitalDelta]t},t1=AbsoluteTime[];\[CapitalDelta]t= UnitConvert[ Quantity[N[t1-t4v], "Seconds" ], "Hours" ];
Print[ "4 vortices loop timing \[CapitalDelta]t = ",IntegerPart[\[CapitalDelta]t],IntegerPart@UnitConvert[FractionalPart[\[CapitalDelta]t], "Minutes" ]   ] ];
