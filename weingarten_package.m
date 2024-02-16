(* ::Package:: *)

BeginPackage["Weingarten`"];
(*names to exports*)


murnNaka::usage="murnNaka[lambda_List,mu_List]: Returns \!\(\*SubsuperscriptBox[\(\[Chi]\), \(\[Mu]\), \(\[Lambda]\)]\)";

badMapping::usage="badMapping[tableaux_List, mu_List]: Takes in a single tableaux and a character mu, and returns false if it can't be filled appropriately with the elements of Mu.";

onlyBorderstrip::usage="onlyBorderstrip[tableauxList_List, mu_List]: Takes a list of tableaux, which then have the tableaux that are not compatible with mu deleted from them.";

SSYT::usage="SSYT[seedList_List,lambda_List,mu_List]: Takes in a list and builds all eligible semi-standard young tableaux in an overpopulous list." ;

derivativeTableaux::usage="derivativeTableaux[tableaux_List,addUnit_,lambda_List]: Takes a single tableau and adds the selected number to it's contents in a way that keeps it semi-standard. Prints all tableaux that fit this rule.";

getClass::usage = "getClass[p_,cycle_] returns the class in \!\(\*SubscriptBox[\(S\), \(p\)]\) from the integer p and a cycle in \!\(\*SubscriptBox[\(S\), \(p\)]\).";

snDimension::usage = "snDimension[partition_] returns the dimension of the representation partition in Sn";

udDimension::usage = "udDimension[ppartition_,d_] returns the dimension of irrep `ppartition` of U(d).";

gClass::usage = "gclass[partition_] returns the number of elements in the class `partition`.";

eWg::usage = "eWg[sigma_,lengthofproduct_,d_] returns the Weingarten function for the product of `lengthofproduct` U's and \!\(\*SuperscriptBox[\(U\), \(*\)]\)'s of U(d) matrices unchanged under the permutation cycle sigma.";

cWg::usage = "cWg[class_,d_] returns the Weingarten function for a product of U's and \!\(\*SuperscriptBox[\(U\), \(*\)]\)'s of U(d) matrices unchanged under permutations in the class `class`.";





Begin["`private`"];
(*functions*)

murnNaka[lambda_List,mu_List]:=Module[{d1,d2,tableauxList,i,j,k,height,weight,seedList},
(*Semi-lazy check to see if the partitions will even work together. Add validPartition if seriously needed*)
d1=Length[lambda];
d2=Length[mu];
If[Total[lambda]!=Total[mu],Return[0]];

seedList={ConstantArray[0,{d1,lambda[[1]]}]};(*Creates the list used to seed the SSYT program*)
tableauxList=SSYT[seedList,lambda,mu];(*Generates an overpopulous list of SSYT*)
tableauxList=onlyBorderstrip[tableauxList,mu];(*Removes all tableaux that aren't borderstrip tableaux*)

For[i=1,i<=Length[tableauxList],i++, (*Removes all tableaux that have 2x2 squares, this should probably be modified so it is included in onlyBorderstrip*)
For[j=1,j<d1,j++,
For[k=1,k<lambda[[1]],k++,
If[i==0,Break[]];(*If future problems are found, blame this first.*)
If[(tableauxList[[i,j,k]]!=0)&&(tableauxList[[i,j,k]]===tableauxList[[i,j+1,k]])&&(tableauxList[[i,j,k]]===tableauxList[[i,j,k+1]])&&(tableauxList[[i,j,k]]==tableauxList[[i,j+1,k+1]]),
tableauxList=Delete[tableauxList,i];
i--;
];
];
];
];

(*Counts the heights and calculates the weight*)
For[i=1,i<=Length[tableauxList],i++,
For[j=1,j<=d2,j++,
height[i,j]=-1;
For[k=d1,k>=1,k--,
If[Count[tableauxList[[i,k]],j]>0,
height[i,j]++;
];
];
];
weight[i]=(-1)^Sum[height[i,jj],{jj,1,d2}];
];

(*Sums the weights for the character*)
Sum[weight[ii],{ii,1,Length[tableauxList]}]
]

(*
badMapping
Dependencies: None
Known Issues: As detailed in murnNaka
Improvements: As detailed in murnNaka
*)

badMapping[tableaux_List,mu_List]:=Module[{M,i,d,count,badMap=False,j,k,masterd={}},
(*Print[tableaux//MatrixForm];*)
For[i=1,i<=Length[mu],i++,
M[i]=Position[tableaux,i];
count[i]=Count[Flatten[tableaux],i];
If[Length[M[i]]==1,M[i]=Append[M[i],M[i][[1]]+{1,0}];];
M[i]=ArrayReshape[Tuples[M[i],2],{count[i],count[i],2,2}];
For[j=1,j<=count[i],j++,
For[k=1,k<=count[i],k++,
M[i]=ReplacePart[M[i],{j,k}->EuclideanDistance[M[i][[j,k,1]],M[i][[j,k,2]]]];
];
];
d[i]={};

For[j=1,j<=count[i],j++,
d[i]=Append[d[i],Count[M[i][[j]],1]]
                  ];

For[j=1,j<=Length[d[i]],j++,
If[count[i]>= 2&& d[i][[j]]==0,badMap=True];
];

d[i]=Count[d[i],1];
If[d[i]==0&&count[i]!=1,d[i]=3];(*CURRENT ERROR TRAPPING FOR ROGUE 1x1s, DONT TRUST IT*)
masterd=Append[masterd,d[i]];
];
If[Max[masterd]>2,badMap=True];
badMap
]

(*
onlyBorderstrip
Dependencies: badMapping
Known Issues: As in murnNaka
Improvements: As in murnNaka
*)

onlyBorderstrip[tableauxList_List,mu_List]:=Module[{i,copyList},
copyList=tableauxList;
(*For each tableaux in the list, delete the one that has 'bad mappings'*)
For[i=1,i<=Length[copyList],i++,
If[badMapping[copyList[[i]],mu],
copyList=Delete[copyList,i];
i--;
];
];
copyList
]

(*
SSYT, semi-standard young tableaux
Dependencies: derivativeTableaux, SSYT
Known Issues: Recursive
Improvements: Can the recursive definition be removed, or made faster. Replace Append with Reap and Sow.
*)

SSYT[seedList_List,lambda_List,mu_List]:=Module[{return,addUnit,i,newmu,newTableauxList},
If[Total[mu]==0,return=seedList];
If[Total[mu]!=0,
newmu=mu;
newTableauxList={};
For[i=1,i<=Length[mu],i++,If[mu[[i]]!=0,addUnit=i;newmu[[i]]--;Break[];];];
For[i=1,i<=Length[seedList],i++,
newTableauxList=Append[newTableauxList,derivativeTableaux[seedList[[i]],addUnit,lambda]];
];
newTableauxList=Flatten[newTableauxList,1];
newTableauxList=DeleteDuplicates[newTableauxList];
return=SSYT[newTableauxList,lambda,newmu]
];
return
]

(*
derivativeTableaux
Dependencies: None
Known Issues: None
Improvements: Replace Append with Reap and Sow.
*)

derivativeTableaux[tableaux_List,addUnit_,lambda_List]:=Module[{rows,cols,i,j,tableauxDerivatives={}},
rows=Length[lambda];
cols=lambda[[1]];
For[i=1,i<=rows,i++,
For[j=1,j<=lambda[[i]],j++,
If[(tableaux[[i,j]]==0)&&(
(i==1&&j==1&&addUnit==1)
||
(j==1&&tableaux[[i-1,1]]<=addUnit&&tableaux[[i-1,1]]!=0)
||
(i==1&&tableaux[[1,j-1]]<=addUnit&&tableaux[[1,j-1]]!=0)
||
(i!=1&&j!=1&&tableaux[[i-1,j]]<=addUnit&&tableaux[[i,j-1]]<=addUnit&&tableaux[[i-1,j]]!=0&&tableaux[[i,j-1]]!=0)
),
tableauxDerivatives=Append[tableauxDerivatives,ReplacePart[tableaux,{i,j}-> addUnit]];
];
];
];
tableauxDerivatives
]

(*
tableaux, print, tableux, printer, tableauxPrinter
Dependencies: None
Known Issues: None
Improvements: None
*)
(* tableauxPrinter::usage="tableauxPrinter[tableauxList_List]: Prints a list of tableaux. How conveninant." *)
tableauxPrinter[tableauxList_List]:=Module[{i},
For[i=1,i<=Length[tableauxList],i++,
Print[tableauxList[[i]]//MatrixForm];
];
]

(* d as in U(d)
p is the length of the product *)

getClass::orderTooLow="The order you have provided (`1`) is too low. It must be an integer greater than 2.";
getClass::badArgs="Bad arguments. getClass takes exactly 2 arguments: an 'order' and a 'cycle'. The order must be an integer greater than 2 and the cycle must be a valid permutation cycle.";

getClass[p_?IntegerQ,cycle_?PermutationCyclesQ]:= Module[{listofcycles,shortlistofcycles,padding},
listofcycles=Flatten[List@@cycle,1];
shortlistofcycles=Table[Length[listofcycles[[k]]],{k,1,Length[listofcycles]}];
padding=p-Sum[shortlistofcycles[[k]],{k,1,Length[shortlistofcycles]}];
If[padding!=0,AppendTo[shortlistofcycles,Table[1,{q,1,padding}]]];
Reverse[Sort[Flatten[shortlistofcycles]]]]

(* 
character::usage = "character[partition_, class_] returns the character for irrep partition of any element in class";

character[partition_,class_]:= murnNaka[partition,class];  *)


snDimension[partition_]:=Module[{numerator,kk,denominator,enn},
enn=Sum[partition[[k]],{k,1,Length[partition]}];
kk=Length[partition];
numerator=Product[(partition[[i]]-partition[[j]]+j-i),{i,1,kk},{j,i+1,kk}];
denominator=Product[(partition[[i]]+kk-i)!,{i,1,kk}];
enn! numerator/denominator
]


 udDimension[ppartition_,d_]:= (* get dimension of irrep ppartition of U(d) as a function of d *)
Module[{listofparts,numerator,rpartition,denominator,conjugatepartition,hooks},  
listofparts=Table[Product[(d+k-i),{k,1,ppartition[[i]]}],{i,1,Length[ppartition]}];
numerator=Times@@listofparts;
(* hook lengths *)
rpartition=Reverse[ppartition];
conjugatepartition=Table[Sum[If[i<=rpartition[[j]],1,0],{j,1,Length[ppartition]}],{i,1,ppartition[[1]]}];
hooks=Flatten[Table[1+ppartition[[i]]-j+conjugatepartition[[j]]-i,{i,1,Length[ppartition]},{j,1,ppartition[[i]]}],1];
denominator=Times@@hooks;
numerator/denominator
]

(* ---  *)

eWg[sigma_?PermutationCyclesQ,lengthofproduct_,d_]:= Module[{class,dimensions,partitionlist},
class=getClass[lengthofproduct,sigma];
partitionlist=IntegerPartitions[lengthofproduct];
dimensions=Table[snDimension[partitionlist[[k]]],{k,1,Length[partitionlist]}];
Sum[dimensions[[k]]^2 murnNaka[partitionlist[[k]],class]/udDimension[partitionlist[[k]],d]/((lengthofproduct!)^2),{k,1,Length[partitionlist]}]
]

cWg[class_,d_]:= Module[{dimensions,lengthofproduct,partitionlist},
lengthofproduct=Total[class];
partitionlist=IntegerPartitions[lengthofproduct];
dimensions=Table[snDimension[partitionlist[[k]]],{k,1,Length[partitionlist]}];
Sum[dimensions[[k]]^2 murnNaka[partitionlist[[k]],class]/udDimension[partitionlist[[k]],d]/((lengthofproduct!)^2),{k,1,Length[partitionlist]}]//Together
]

(* ---  *)

gClass[partition_]:=Module[{enn,nparts,tableofcycles},  (* number of elements in a class *)
enn=Total[partition];
nparts=Length[partition];
tableofcycles=Tally[partition];
enn!/Product[tableofcycles[[k,1]]^tableofcycles[[k,2]],{k,1,Length[tableofcycles]}]/Product[tableofcycles[[k,2]]!,{k,1,Length[tableofcycles]}]]

  
End[];
EndPackage[]
