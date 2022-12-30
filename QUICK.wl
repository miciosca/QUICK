(* ::Package:: *)

If[MemberQ[$Packages,"QUICK`"],Print["The package QUICK is already loaded!"];Abort[]]


BeginPackage["QUICK`"]


Print["-*-*-*-*-",Style[#,Bold]&@"QUICK","-1.0-*-*-*-*-\n",Style[#,Bold]&@"
QUI","ver ",Style[#,Bold]&@"C","orrelators ",Style[#,Bold]&@"K","it.\n
Author: Michelangelo Preti, King's College London\n
-*-*-*-*-*-*-*-*-*-*-*-*-*-\n
QUICK loaded! The new provided commands are:\n 
- ComputeOO, ComputeW and ComputeWO\n 
- NormalOrderedOP, SphereCorrelatorOO, SphereCorrelatorWO, GramSchmidtCoeff\n 
- QUICKsaveData, QUICKsaveResults"]


ComputeOO::usage="ComputeOO[Quiver,PerturbativeOrder][{NodeO1,dimO1},{NodeO2,dimO2}]

It computes the correlator \!\(\*SuperscriptBox[SubscriptBox[\(G\), \(dimO1\)], \((NodeO1, NodeO2)\)]\)(\!\(\*SubscriptBox[\(\[Lambda]\), \(1\)]\),...,\!\(\*SubscriptBox[\(\[Lambda]\), \(Quiver\)]\),N) = \[LeftAngleBracket]:\!\(\*SuperscriptBox[SubscriptBox[\(O\), \(dimO1\)], \((NodeO1)\)]\): :\!\(\*SuperscriptBox[SubscriptBox[\(O\), \(dimO2\)], \((NodeO2)\)]\):\!\(\*SubscriptBox[\(\[RightAngleBracket]\), \(Quiver\)]\) up to order \!\(\*SuperscriptBox[\(\[Lambda]\), \(PerturbativeOrder/2\)]\).
dimO1 and dimO2 are lists, the other inputs are non-negative integers.
In this case dimO2=dimO1 since normal-ordered operators are orthogonal."
ComputeW::usage="ComputeW[Quiver,PerturbativeOrder][NodesWL]
 
It computes the expectation value of Wilson loops \!\(\*SuperscriptBox[SubscriptBox[\(w\), \(NodesWL\)], \((Quiver)\)]\)(\!\(\*SubscriptBox[\(\[Lambda]\), \(1\)]\), ..., \!\(\*SubscriptBox[\(\[Lambda]\), \(Quiver\)]\),N)=\[LeftAngleBracket]\!\(\*SubscriptBox[\(W\), \(NodesWL\)]\)\!\(\*SubscriptBox[\(\[RightAngleBracket]\), \(Quiver\)]\) up to order \!\(\*SuperscriptBox[\(\[Lambda]\), \(PerturbativeOrder/2\)]\).
NodesWL is a list, the other inputs are non-negative integers."
ComputeWO::usage="ComputeW[Quiver,PerturbativeOrder][NodesWL,{NodeO,dimO}]

It computes the correlator \!\(\*SuperscriptBox[SubscriptBox[\(A\), \(dimO\)], \((NodesWL, NodeO)\)]\)(\!\(\*SubscriptBox[\(\[Lambda]\), \(1\)]\), ..., \!\(\*SubscriptBox[\(\[Lambda]\), \(Quiver\)]\), N )=\[LeftAngleBracket]\!\(\*SubscriptBox[\(W\), \(NodesWL\)]\) :\!\(\*SuperscriptBox[SubscriptBox[\(O\), \(dimO\)], \((NodeO)\)]\):\!\(\*SubscriptBox[\(\[RightAngleBracket]\), \(Quiver\)]\) up to order \!\(\*SuperscriptBox[\(\[Lambda]\), \(PerturbativeOrder/2\)]\).
dimO and NodesWL are lists, the other inputs are non-negative integers."
NormalOrderedOP::usage="NormalOrderedOP[Quiver][Node,dim]

It gives the normal ordered operator :\!\(\*SuperscriptBox[SubscriptBox[\(O\), \(dim\)], \((Node)\)]\): in terms of the basis of operators \!\(\*OverscriptBox[\(O\), \(~\)]\) by means of (2.34) in the theory \!\(\*SubscriptBox[\(A\), \(Quiver - 1\)]\).
dim is a list, the other inputs are non-negative integers."
SphereCorrelatorOO::usage="SphereCorrelatorOO[Quiver,PerturbativeOrder][{NodeO1,dimO1},{NodeO2,dimO2}]

It computes the multi matrix model correlator \[LeftAngleBracket]\!\(\*SuperscriptBox[SubscriptBox[OverscriptBox[\(O\), \(~\)], \(dimO1\)], \((NodeO1)\)]\)  \!\(\*SuperscriptBox[SubscriptBox[OverscriptBox[\(O\), \(~\)], \(dimO2\)], \((NodeO2)\)]\)\!\(\*SubscriptBox[\(\[RightAngleBracket]\), \(Quiver\)]\) up to \!\(\*SuperscriptBox[\(\[Lambda]\), \(PerturbativeOrder/2\)]\).
dimO1 and dimO2 are lists, the other inputs are non-negative integers."
SphereCorrelatorWO::usage="SphereCorrelatorWO[Quiver,PerturbativeOrder][NodesWL,{NodeO,dimO}]

It computes the multi matrix model correlator \[LeftAngleBracket]\!\(\*SubscriptBox[\(W\), \(NodesWL\)]\) \!\(\*SuperscriptBox[SubscriptBox[OverscriptBox[\(O\), \(~\)], \(dimO\)], \((NodeO)\)]\)\!\(\*SubscriptBox[\(\[RightAngleBracket]\), \(Quiver\)]\) up to \!\(\*SuperscriptBox[\(\[Lambda]\), \(PerturbativeOrder/2\)]\).
dimO and NodesWL are lists, the other inputs are non-negative integers."
GramSchmidtCoeff::usage="GramSchmidtCoeff[Quiver,PerturbativeOrder][Node,dim]

It computes the Gram-Schmidt coefficients appearing in the expansion of the operator :\!\(\*SuperscriptBox[SubscriptBox[\(O\), \(dim\)], \((Node)\)]\): in the theory \!\(\*SubscriptBox[\(A\), \(Quiver \[Minus] 1\)]\) up to \!\(\*SuperscriptBox[\(\[Lambda]\), \(PerturbativeOrder/2\)]\).
dim is a list, the other inputs are non-negative integers."
QUICKsaveData::usage="QUICKsaveData[\"filename.mx\"]

It generates a file filename.mx containing the list of the internal functions already computed in the current Mathematica session.
This file can be loaded and then updated at any new session to speed-up the computational time."
QUICKsaveResults::usage="QUICKsaveResults[\"filename.mx\"]

It generates a file filename.mx containing the list of all the perturbative expansions computed in the current Mathematica session.
This file can be loaded and then updated at any new session to speed-up the computational time."


Begin["`Private`"]


(* ::Subsection:: *)
(*Notation & Initial conditions*)


$PrePrint=TraditionalForm;


Global`$OrbifoldPoint=False;
Global`$LargeN=False;
Global`$UNgroup=False;
Global`$TranscendentalExp=True;


MakeBoxes[Global`\[Zeta][b__],TraditionalForm]:=SubscriptBox[StyleBox["\[Zeta]",FontColor->Red],MakeBoxes[Grid[{{b}}],TraditionalForm]];


MakeBoxes[Global`\[Lambda][b__],TraditionalForm]:=SubscriptBox[StyleBox["\[Lambda]",FontColor->Blue],MakeBoxes[Grid[{{b}}],TraditionalForm]];


SetAttributes[Global`tt,Orderless]
SetAttributes[Global`ttu,Orderless]
SetAttributes[c,Orderless]


(* ::Subsection:: *)
(*Functions*)


Action[q_,ncut_]:=If[q===1,Sum[-((-1)^(n+k)/n) (Global`\[Lambda][1]/(8\[Pi]^2 Global`\[CapitalNu]))^n Global`\[Zeta][2n-1]Binomial[2n,k]tr[X[1]^(2n-k)]tr[X[1]^k],{n,2,ncut},{k,2,2n-2}]/.If[Global`$UNgroup,{},tr[X[p_]]->0],Sum[-((-1)^(n+k)/n) (1/(8\[Pi]^2))^n Global`\[Zeta][2n-1]Binomial[2n,k]Sum[(Global`\[Lambda][i]^n/Global`\[CapitalNu]^(n) tr[X[i]^(2n-k)]tr[X[i]^k]-(Global`\[Lambda][i]^(n-k/2) Global`\[Lambda][i+1]^(k/2))/Global`\[CapitalNu]^(n) tr[X[i]^(2n-k)]tr[X[i+1]^k]),{i,1,q}],{n,2,ncut},{k,0,2n}]/.If[Global`$UNgroup,{},tr[X[p_]]->0]/.X[q+1]->X[1]/.If[Global`$OrbifoldPoint,Global`\[Lambda][l_]:>Global`\[Lambda][1],Global`\[Lambda][q+1]->Global`\[Lambda][1]]]


TwoPF[q_,ord_][{a1_,a2_},{b1_,b2_}]:=Module[{action2pf=Module[{actionexp=Series[Exp[If[EvenQ[Plus@@a2]&&EvenQ[Plus@@b2]&&Global`$LargeN,Action[q,ord/2]/.tr[X[a_]^b_.]/;OddQ[b]:>0,Action[q,ord/2]]]/.Global`\[Lambda][l_]:>Global`\[Epsilon]^2 Global`\[Lambda][l],{Global`\[Epsilon],0,ord},Assumptions->{Global`\[Epsilon]>0,Global`\[Lambda][__]>0}]//PowerExpand},If[Head[actionexp]===SeriesData,actionexp,actionexp+SeriesData[Global`\[Epsilon],0,{},0,ord+1,1]]]},((1/((action2pf//ExpandAll)/.tr[X[p_]^a_.]^b_.:>(\[Tau][p]@@Table[a,{i,1,b}])//.\[Tau][p_][a___]\[Tau][p_][b___]:>\[Tau][p]@@Sort[{a,b}]))(((Times@@(tr[X[a1]^#]&/@a2))(Times@@(tr[X[b1]^#]&/@b2))action2pf//ExpandAll)/.tr[X[p_]^a_.]^b_.:>(\[Tau][p]@@Table[a,{i,1,b}])//.\[Tau][p_][a___]\[Tau][p_][b___]:>\[Tau][p]@@Sort[{a,b}]))/.\[Tau][p_][a___]/;OddQ[Total[List@a]]:>0/.If[Global`$UNgroup,{},\[Tau][p_][a___]/;!FreeQ[List@a,1]:>0]]


OnePF[q_,ord_][{a1_,a2_}]:=Module[{action1pf=Module[{actionexp=Series[Exp[If[EvenQ[Plus@@a2]&&Global`$LargeN,Action[q,ord/2]/.tr[X[a_]^b_.]/;OddQ[b]:>0,Action[q,ord/2]]]/.Global`\[Lambda][l_]:>Global`\[Epsilon]^2 Global`\[Lambda][l],{Global`\[Epsilon],0,ord},Assumptions->{Global`\[Epsilon]>0,Global`\[Lambda][__]>0}]//PowerExpand},If[Head[actionexp]===SeriesData,actionexp,actionexp+SeriesData[Global`\[Epsilon],0,{},0,ord+1,1]]]},(1/((action1pf//ExpandAll)/.tr[X[p_]^a_.]^b_.:>(\[Tau][p]@@Table[a,{i,1,b}])//.\[Tau][p_][a___]\[Tau][p_][b___]:>\[Tau][p]@@Sort[{a,b}]) (((Times@@(tr[X[a1]^#]&/@a2))action1pf//ExpandAll)/.tr[X[p_]^a_.]^b_.:>(\[Tau][p]@@Table[a,{i,1,b}])//.\[Tau][p_][a___]\[Tau][p_][b___]:>\[Tau][p]@@Sort[{a,b}]))/.\[Tau][p_][a___]/;OddQ[Total[List@a]]:>0/.If[Global`$UNgroup,{},\[Tau][p_][a___]/;!FreeQ[List@a,1]:>0]]


NormalOrdered2PF[Quiver_,ord_][{a1_,a2_},{b1_,b2_}]:=If[EvenQ[Plus@@a2]&&EvenQ[Plus@@b2],TwoPF[Quiver,ord][{a1,a2},{b1,b2}]-OnePF[Quiver,ord][{a1,a2}]OnePF[Quiver,ord][{b1,b2}],TwoPF[Quiver,ord][{a1,a2},{b1,b2}]]


\[Tau]solve[A_]:=A/.\[Tau][p_][a___]:>ttsolve[Global`tt[a]]


ttsolve[A_]:=A/.Global`tt[]->(Global`tt[]=1)/.Global`tt[a___]/;OddQ[Total[List@a]]:>0//.Global`tt[a___]/;MemberQ[List@a,1]:>0/.Global`tt[a___]/;MemberQ[List@a,0]:>(Global`\[CapitalNu]^Count[List@a,0] ttsolve[Global`tt@@DeleteCases[List@a,0]]//Simplify)/.Global`tt[a___]:>If[!Or@@(NumberQ[#]&/@{a}),Global`TT[a],(Global`tt[a]=ttsolve[(Module[{list=List@a},1/2 Sum[Global`tt@@MapAt[#-m-2&,Prepend[list,m],2]-1/Global`\[CapitalNu] Global`tt@@MapAt[#-2&,list,1],{m,0,list[[1]]-2}]+Sum[list[[1+m]]/2 (Global`tt@@MapAt[#+list[[1]]-2&,Drop[list,1],m]-1/Global`\[CapitalNu] Global`tt@@Insert[MapAt[#-1&,Drop[list,1],m],list[[1]]-1,m]),{m,1,Length[list]-1}]]//Simplify)])]//Simplify


\[Tau]usolve[A_]:=A/.\[Tau][p_][a___]:>ttusolve[Global`ttu[a]]


ttusolve[A_]:=A/.Global`ttu[]->(Global`ttu[]=1)/.Global`ttu[a___]/;OddQ[Total[List@a]]:>0/.Global`ttu[a___]/;MemberQ[List@a,0]:>(Global`\[CapitalNu]^Count[List@a,0] ttusolve[Global`ttu@@DeleteCases[List@a,0]]//Simplify)/.Global`ttu[a___]:>If[!Or@@(NumberQ[#]&/@{a}),Global`TT[a],(Global`ttu[a]=ttusolve[(Module[{list=List@a},1/2 Sum[Global`ttu@@MapAt[#-m-2&,Prepend[list,m],2](*-(1/Global`\[CapitalNu])Global`tt@@MapAt[#-2&,list,1]*),{m,0,list[[1]]-2}]+Sum[list[[1+m]]/2 (Global`ttu@@MapAt[#+list[[1]]-2&,Drop[list,1],m](*-(1/Global`\[CapitalNu])Global`tt@@Insert[MapAt[#-1&,Drop[list,1],m],list[[1]]-1,m]*)),{m,1,Length[list]-1}]]//Simplify)])]//Simplify


\[Sigma]solve[A_]:=A/.s[node_][nn_][B___]:>sssolve[node][Global`ss[nn][Sequence@@Sort[{B}]]]


sssolve[ii_][B___]:=B/.Global`ss[nn_][A___]:>(Global`ss[nn][A]=((Plus@@Module[{listcoeff=CoefficientList[(1/(Product[(Global`\[Lambda][1]/(2 Global`\[CapitalNu]))^(Global`j[l]/2) 1/Global`j[l]!,{l,1,nn}]Global`TT[Sequence@@Table[Global`j[l],{l,1,nn}]]) (((Product[(Global`\[Lambda][1]/(2 Global`\[CapitalNu]))^(Global`k[l]/2) 1/Global`k[l]!,{l,1,nn}]If[A===0,\[Tau][1][Sequence@@Table[Global`k[l],{l,1,nn}]],\[Tau][1][A,Sequence@@Table[Global`k[l],{l,1,nn}]]]//\[Tau]solve)//Collect[#,Global`TT[__]]&)//.coeff_ Global`TT[QQ___,Global`k[w_]+a_.,PP___]:>(coeff/.Global`k[w]->Global`j[w]-a)Global`TT[QQ,Global`j[w],PP])//FunctionExpand)/.Global`j[nn]->-Sum[Global`j[l],{l,1,nn-1}]+sum,sum]//Simplify},(Table[ 2^n Der[Global`W[nn],{Global`\[Lambda][1],n}],{n,0,Length[listcoeff]-1}]listcoeff//Expand)//.Global`j[a_]^k_. Der[QQ_,b_]:>Der[Global`j[a]^k QQ,b]/.Der[QQ_,b_]/;!FreeQ[QQ,Global`j]:>Der[Global`\[CapitalSigma][QQ],b]//Simplify]/.Der[Global`W[nn],{Global`\[Lambda][p_],n_}]:>Nest[Simplify[Global`\[Lambda][p] D[#,Global`\[Lambda][p]]]&,Global`W[nn][Global`\[Lambda][p]],n]/.Der[Global`\[CapitalSigma][AAA_],{Global`\[Lambda][p_],n_}]:>Nest[Simplify[Global`\[Lambda][p] D[#,Global`\[Lambda][p]]]&,Global`\[CapitalSigma][AAA][Global`\[Lambda][p]],n])//Collect[#,Derivative[___][Global`W[M_]][Global`\[Lambda][MM_]],Simplify]&))/.If[Global`$OrbifoldPoint,{},Global`\[Lambda][p_]:>Global`\[Lambda][ii]]


\[Sigma]usolve[A_]:=A/.s[node_][nn_][B___]:>ssusolve[node][Global`ssu[nn][Sequence@@Sort[{B}]]]


ssusolve[ii_][B___]:=B/.Global`ssu[nn_][A___]:>(Global`ssu[nn][A]=((Plus@@Module[{listcoeff=CoefficientList[(1/(Product[(Global`\[Lambda][1]/(2 Global`\[CapitalNu]))^(Global`j[l]/2) 1/Global`j[l]!,{l,1,nn}]Global`TT[Sequence@@Table[Global`j[l],{l,1,nn}]]) (((Product[(Global`\[Lambda][1]/(2 Global`\[CapitalNu]))^(Global`k[l]/2) 1/Global`k[l]!,{l,1,nn}]If[A===0,\[Tau][1][Sequence@@Table[Global`k[l],{l,1,nn}]],\[Tau][1][A,Sequence@@Table[Global`k[l],{l,1,nn}]]]//\[Tau]usolve)//Collect[#,Global`TT[__]]&)//.coeff_ Global`TT[QQ___,Global`k[w_]+a_.,PP___]:>(coeff/.Global`k[w]->Global`j[w]-a)Global`TT[QQ,Global`j[w],PP])//FunctionExpand)/.Global`j[nn]->-Sum[Global`j[l],{l,1,nn-1}]+sum,sum]//Simplify},(Table[ 2^n Der[Global`W[nn],{Global`\[Lambda][1],n}],{n,0,Length[listcoeff]-1}]listcoeff//Expand)//.Global`j[a_]^k_. Der[QQ_,b_]:>Der[Global`j[a]^k QQ,b]/.Der[QQ_,b_]/;!FreeQ[QQ,Global`j]:>Der[Global`\[CapitalSigma][QQ],b]//Simplify]/.Der[Global`W[nn],{Global`\[Lambda][p_],n_}]:>Nest[Simplify[Global`\[Lambda][p] D[#,Global`\[Lambda][p]]]&,Global`W[nn][Global`\[Lambda][p]],n]/.Der[Global`\[CapitalSigma][AAA_],{Global`\[Lambda][p_],n_}]:>Nest[Simplify[Global`\[Lambda][p] D[#,Global`\[Lambda][p]]]&,Global`\[CapitalSigma][AAA][Global`\[Lambda][p]],n])//Collect[#,Derivative[___][Global`W[M_]][Global`\[Lambda][MM_]],Simplify]&))/.If[Global`$OrbifoldPoint,{},Global`\[Lambda][p_]:>Global`\[Lambda][ii]]


distance[Quiver_,I_,J_]:=Abs[Abs[Abs[J-I]-Quiver/2]-Quiver/2]


which\[Alpha][qq_]:=(Cases[qq,If[Global`$OrbifoldPoint,Global`\[Alpha][_,_][_,_],Global`\[Alpha][_,_,_][_,_]],\[Infinity]]//Union);


whichc[qq_]:=(Cases[{qq},c[_,_],\[Infinity]]//Union);


partitions[n_,k_]:=(Composition[Union,Permutations,PadRight[#,k]&]/@IntegerPartitions[n,k]//Flatten[#,1]&//Sort)


deletemultitrace[list_,length_]:=Pick[list,UnitStep[Length/@list-(1+length)],0]


listoperators[A___]:=Module[{list=Insert[(Sort/@(Flatten[Table[partitions[Plus[A]-2k,i],{k,1,Plus[A]/2},{i,1,Plus[A]-2k}],2]//If[Global`$UNgroup,DeleteCases[#,0,Infinity],DeleteCases[DeleteCases[#,0,Infinity],{___, 1,___}]]&))//DeleteDuplicates,List[A],1]//SortBy[#,Total]&//Reverse},If[Global`$LargeN,deletemultitrace[list,Length[{A}]],list]]


(*listoperatorsW[A___]:=Insert[(Sort/@(Flatten[Table[partitions[Plus[A]-2k,i],{k,1,Plus[A]/2},{i,1,Plus[A]-2k}],2]//If[Global`$UNgroup,DeleteCases[#,0,Infinity],DeleteCases[DeleteCases[#,0,Infinity],{___, 1,___}]]&))//DeleteDuplicates,List[A],1]//SortBy[#,Total]&//Reverse*)


XX[Quiver_,Node_,dim_]:=Module[{list=listoperators[Sequence@@dim]},If[Global`$OrbifoldPoint,Insert[Table[Global`\[Alpha][Quiver,distance[Quiver,Node,s]][dim,list[[k]]]\[Chi][s,list[[k]]],{k,2,Length[list]},{s,1,Quiver}]//Flatten,\[Chi][Node,list[[1]]],1],Insert[Table[Global`\[Alpha][Quiver,Node,s][dim,list[[k]]]\[Chi][s,list[[k]]],{k,2,Length[list]},{s,1,Quiver}]//Flatten,\[Chi][Node,list[[1]]],1]]]


(*XXW[Quiver_,Node_,dim_]:=Module[{list=listoperatorsW[Sequence@@dim]},If[Global`$OrbifoldPoint,Insert[Table[Global`\[Alpha][Quiver,distance[Quiver,Node,s]][dim,list[[k]]]\[Chi][s,list[[k]]],{k,2,Length[list]},{s,1,Quiver}]//Flatten,\[Chi][Node,list[[1]]],1],Insert[Table[Global`\[Alpha][Quiver,Node,s][dim,list[[k]]]\[Chi][s,list[[k]]],{k,2,Length[list]},{s,1,Quiver}]//Flatten,\[Chi][Node,list[[1]]],1]]]*)


LinSol[Eqs_,Vars_]:=Module[{b,A,Sol,SolRule},{b,A}=CoefficientArrays[Eqs,Vars];
Sol=LinearSolve[A,-b,Method->"CofactorExpansion"];
SolRule={};
Do[SolRule=Union[SolRule,{Vars[[i]]->Sol[[i]]}],{i,1,Length[Vars]}];
SolRule]


GSeqs[Quiver_,Node_,dim_]:=Module[{NOop=XX[Quiver,Node,dim]},Module[{ops=Drop[NOop,1]/.Global`\[Alpha][__][__]:>1},Table[Plus@@NOop/.\[Chi][p1_,d_]:>c[\[Chi][p1,d],ops[[k]]],{k,1,Length[ops]}]]]


(*GSeqsW[Quiver_,Node_,dim_]:=Module[{NOop=XXW[Quiver,Node,dim]},Module[{ops=Drop[NOop,1]/.Global`\[Alpha][__][__]:>1},Table[Plus@@NOop/.\[Chi][p1_,d_]:>c[\[Chi][p1,d],ops[[k]]],{k,1,Length[ops]}]]]*)


seriesOrder[A_]:=A/.SeriesData[a_,b_,c_,d_,order_,e_]:>(order-1)//Quiet


whichcw[qq_]:=(Cases[qq,cw[_,\[Chi][_,_]],\[Infinity]]//Union);


MissingOrders[Quiver_,ord_,j1_]:=With[{orb=Global`$OrbifoldPoint,largeN=Global`$LargeN,un=Global`$UNgroup,texp=Global`$TranscendentalExp},#/.a_/;Negative[#]:>0&/@(ord-(j1/.P[dist_][dim1_,dim2_]:>If[If[$VersionNumber<12.2,!ValueQ[Global`CorrSphere[orb,largeN,un][Quiver,dist][dim1,dim2]],!ValueQ[Global`CorrSphere[orb,largeN,un][Quiver,dist][dim1,dim2],Method->"Legacy"]],-1,seriesOrder[Global`CorrSphere[orb,largeN,un][Quiver,dist][dim1,dim2]]]/.cw[posWL_,\[Chi][aa_,bb_]]:>If[If[$VersionNumber<12.2,!ValueQ[Global`CorrSphereW[orb,largeN,un,texp][Quiver][posWL,{aa,bb}]],!ValueQ[Global`CorrSphereW[orb,largeN,un,texp][Quiver][posWL,{aa,bb}],Method->"Legacy"]],-1,seriesOrder[Global`CorrSphereW[orb,largeN,un,texp][Quiver][posWL,{aa,bb}]]]))]


ParallelCC[Quiver_,ord_][cfromlist_,missingorder_]:=cfromlist/.P[dist_][dim1_,dim2_]:>If[missingorder==0,0,(NormalOrdered2PF[Quiver,ord][{1,dim1},{1+dist,dim2}]/.SeriesData[x_,x0_,coef_,nmin_,nmax_,den_]:>SeriesData[x,x0,Drop[Flatten[{Table[0,{v,1,nmin}],coef,Table[0,{vv,1,nmax-nmin-Length[coef]}]}],nmax-missingorder],nmax-missingorder,nmax,den]//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve])/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,(c//Collect[#,Global`\[Zeta][__],Simplify]&),nmin,nmax,e]//Quiet]/.cw[pWL_,\[Chi][posO_,dim_]]:>If[missingorder==0,0,If[Global`$TranscendentalExp,(NormalOrdered1PFW[Quiver,ord][pWL,{posO,dim}]/.SeriesData[x_,x0_,coef_,nmin_,nmax_,den_]:>SeriesData[x,x0,Drop[Flatten[{Table[0,{v,1,nmin}],coef,Table[0,{vv,1,nmax-nmin-Length[coef]}]}],nmax-missingorder],nmax-missingorder,nmax,den]//If[Global`$UNgroup,\[Sigma]usolve,\[Sigma]solve]//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve]),((NormalOrdered1PFW[Quiver,ord][pWL,{posO,dim}]//\[Sigma]solve\[Lambda]exp[#,ord]&)/.SeriesData[x_,x0_,coef_,nmin_,nmax_,den_]:>SeriesData[x,x0,Drop[Flatten[{Table[0,{v,1,nmin}],coef,Table[0,{vv,1,nmax-nmin-Length[coef]}]}],nmax-missingorder],nmax-missingorder,nmax,den]//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve])]/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,c//If[Global`$TranscendentalExp,Collect[#,Derivative[___][Global`W[A_]][Global`\[Lambda][B_]],Simplify],Collect[#,Global`\[Zeta][__],Simplify]]&,nmin,nmax,e]//Quiet]


ComputeCC[Quiver_,ord_][jj_]:=Module[{j0=Flatten[{whichc[jj],whichcw[jj]}]},Module[{j1=j0/.c[\[Chi][a1_,b1_],\[Chi][a2_,b2_]]:>P[distance[Quiver,a1,a2]][Sequence@@Sort[{b1,b2}]]//DeleteDuplicates},Module[{j2=j1//MissingOrders[Quiver,ord,#]&},Module[{Global`clist=Delete[j1,Position[j2,0]],Global`missing=DeleteCases[j2,0]},Module[{Global`nP=Count[Global`clist,P[_][_,_]],Global`ncw=Count[Global`clist,cw[_,_]]},If[Global`nP===0&&Global`ncw===0,PrintTemporary["All correlators on the sphere at order "<> ToString[ord]<>" are already computed!"];,
print0=PrintTemporary["Number of correlators on the sphere reduced from "<>ToString[Length[j0]]<> " to "<>ToString[Length[j1]]<> " using the symmetries of the quiver."];
If[Length[j1]===Length[Global`clist],print1=PrintTemporary["Parallel computing of "<>ToString[Length[j1]]<>" correlators on the sphere up to order "<>ToString[ord]];,
print1=PrintTemporary["Parallel computing of "<>ToString[Length[Global`clist]]<>" correlators on the sphere up to order "<>ToString[ord]<>"
The remaining "<>ToString[Length[j1]-Length[Global`clist]]<> " are already computed"];];
Module[{indices=Global`clist/.P[dist_][dim1_,dim2_]:>{dist,dim1,dim2}/.cw[a_,\[Chi][b_,c_]]:>{a,b,c},parallel=Module[{Global`progress=0,Global`progressW=0},SetSharedVariable[Global`progress];
SetSharedVariable[Global`progressW];
Monitor[ParallelTable[status[$KernelID]={$KernelID,Global`k,DateString[]};
out=ParallelCC[Quiver,ord][Global`clist[[Global`k]],Global`missing[[Global`k]]];
status[$KernelID]={$KernelID,"DONE!",DateString[]};
If[Global`k<=Global`nP,Global`progress++,Global`progressW++];
out,{Global`k,1,Global`nP+Global`ncw}],Row[{status/@kernels//TableForm,Which[Global`ncw===0,Column[{ToString@Global`progress<>" of "<>ToString[Global`nP],ProgressIndicator[Global`progress,{0,Global`nP}],"\[LeftAngleBracket] \!\(\*SuperscriptBox[SubscriptBox[\(O\), OverscriptBox[\(n\), \(\[Rule]\)]], \(I\)]\) \!\(\*SuperscriptBox[SubscriptBox[\(O\), OverscriptBox[\(p\), \(\[Rule]\)]], \(J\)]\) \[RightAngleBracket]"},Alignment->Center],Global`nP===0,Column[{ToString@Global`progressW<>" of "<>ToString[Global`ncw],ProgressIndicator[Global`progressW,{0,Global`ncw}],"\[LeftAngleBracket] \!\(\*SubscriptBox[\(W\), OverscriptBox[\(p\), \(\[Rule]\)]]\) \!\(\*SuperscriptBox[SubscriptBox[\(O\), OverscriptBox[\(n\), \(\[Rule]\)]], \(I\)]\) \[RightAngleBracket]"},Alignment->Center],True,Sequence@@{Column[{ToString@Global`progress<>" of "<>ToString[Global`nP],ProgressIndicator[Global`progress,{0,Global`nP}],"\[LeftAngleBracket] \!\(\*SuperscriptBox[SubscriptBox[\(O\), OverscriptBox[\(n\), \(\[Rule]\)]], \(I\)]\) \!\(\*SuperscriptBox[SubscriptBox[\(O\), OverscriptBox[\(p\), \(\[Rule]\)]], \(J\)]\) \[RightAngleBracket]"},Alignment->Center],Column[{ToString@Global`progressW<>" of "<>ToString[Global`ncw],ProgressIndicator[Global`progressW,{0,Global`ncw}],"\[LeftAngleBracket] \!\(\*SubscriptBox[\(W\), OverscriptBox[\(p\), \(\[Rule]\)]]\) \!\(\*SuperscriptBox[SubscriptBox[\(O\), OverscriptBox[\(n\), \(\[Rule]\)]], \(I\)]\) \[RightAngleBracket]"},Alignment->Center]}]},"   "]]]},Table[Global`CorrSphere[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][Quiver,indices[[s,1]]][indices[[s,2]],indices[[s,3]]]=Normal[Global`clist[[s]]/.P[dist_][dim1_,dim2_]:>If[Global`missing[[s]]===ord+1,0,Global`CorrSphere[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][Quiver,dist][dim1,dim2]]]+parallel[[s]]//PowerExpand,{s,1,Global`nP}];
Table[Global`CorrSphereW[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup,Global`$TranscendentalExp][Quiver][indices[[s,1]],{indices[[s,2]],indices[[s,3]]}]=Normal[Global`clist[[s]]/.cw[pWL_,\[Chi][posO_,dim_]]:>If[Global`missing[[s]]===ord+1,0,Global`CorrSphereW[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup,Global`$TranscendentalExp][Quiver][pWL,{posO,dim}]]]+parallel[[s]]//PowerExpand,{s,Global`nP+1,Global`nP+Global`ncw}];];
NotebookDelete[{print0,print1}];
PrintTemporary["Correlators on the sphere computed at order "<>ToString[ord]];]]]]]];


crules[Quiver_,AA_]:=AA/.If[Global`$OrbifoldPoint,c[\[Chi][a1_,b1_],\[Chi][a2_,b2_]]:>Global`CorrSphere[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][Quiver,distance[Quiver,a1,a2]][Sequence@@Sort[{b1,b2}]],c[\[Chi][a1_,b1_],\[Chi][a2_,b2_]]:>If[With[{dd=distance[Quiver,a1,a2],orb=Global`$OrbifoldPoint,largeN=Global`$LargeN,un=Global`$UNgroup},If[$VersionNumber<12.2,ValueQ[Global`CorrSphere[orb,largeN,un][Quiver,dd][b1,b2]],ValueQ[Global`CorrSphere[orb,largeN,un][Quiver,dd][b1,b2],Method->"Legacy"]]],If[a2-a1<=Floor[Quiver/2],Global`CorrSphere[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][Quiver,distance[Quiver,a1,a2]][b1,b2]/.(Table[Global`\[Lambda][l]->Global`\[Lambda][l+a1-1],{l,1,Quiver}]/.Global`\[Lambda][aa_]/;aa>Quiver:>Global`\[Lambda][aa-Quiver]),Global`CorrSphere[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][Quiver,distance[Quiver,a1,a2]][b1,b2]/.Flatten[Table[Global`\[Lambda][l]->Global`\[Lambda][1+a1-l],{l,1,Quiver}]/.Global`\[Lambda][0]:>Global`\[Lambda][Quiver]/.Global`\[Lambda][aa_]/;Negative[aa]:>Global`\[Lambda][aa+Quiver]/.Rule[Global`\[Lambda][aa_],Global`\[Lambda][aa_]]:>{}]],If[a2-a1<=Floor[Quiver/2],Global`CorrSphere[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][Quiver,distance[Quiver,a2,a1]][b2,b1] /.Flatten[Table[Global`\[Lambda][l]->Global`\[Lambda][1+a2-l],{l,1,Quiver}]/.Global`\[Lambda][0]:>Global`\[Lambda][Quiver]/.Global`\[Lambda][aa_]/;Negative[aa]:>Global`\[Lambda][aa+Quiver]/.Rule[Global`\[Lambda][aa_],Global`\[Lambda][aa_]]:>{}],Global`CorrSphere[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][Quiver,distance[Quiver,a2,a1]][b2,b1]/.(Table[Global`\[Lambda][l]->Global`\[Lambda][l+a2-1],{l,1,Quiver}]/.Global`\[Lambda][aa_]/;aa>Quiver:>Global`\[Lambda][aa-Quiver]) ]]]


OperatorsByDimension[\[CapitalDelta]_]:=(Sort[#,Greater]&/@(Table[(partitions[\[CapitalDelta],k]//DeleteCases[#,{___, 0,___}]&)//DeleteCases[#,{___, 1,___}]&,{k,1,\[CapitalDelta]/2}]//Flatten[#,1]&))//DeleteDuplicates


Ansatz[A_,cut_]:=A/.Global`\[Alpha][a___][b_,c_]:>Sum[Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][a][b,c,2k]Global`\[Epsilon]^(2k),{k,0,cut}]


whichGScoeff[qq_]:=(Cases[qq,If[Global`$OrbifoldPoint,Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][_,_][_,_,_],Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][_,_,_][_,_,_]],\[Infinity]]//Union);


Compute\[Alpha]temp[Quiver_,ord_,start_,eqs1_]:=For[i=start,i<=Floor[ord/2],i++,tempprint=PrintTemporary["Order "<>ToString[2i]<>"     ..."];
Module[{tempeq=Series[Ansatz[eqs1,i]//crules[Quiver,#]&,{Global`\[Epsilon],0,2i}]//Simplify//Normal},Module[{\[Alpha]=whichGScoeff[tempeq],sol=LinSol[tempeq,whichGScoeff[tempeq]]//Simplify},Module[{indices=\[Alpha]/.Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][a___][b___]:>{{a},{b}}},Table[Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][Sequence@@indices[[k,1]]][Sequence@@indices[[k,2]]]=\[Alpha][[k]]/.sol,{k,1,Length[\[Alpha]]}]]];];
NotebookDelete[tempprint];
ToExpression["tempord"<>ToString[2i]<>"="<>ToString[PrintTemporary["Order "<>ToString[2i]<>"     DONE!"]]];]


Compute\[Alpha][Quiver_,ord_,eqs1_]:=If[Length[which\[Alpha][eqs1]]===0,PrintTemporary["Gram-Schmidt procedure is not needed."];,Module[{prevord=Table[which\[Alpha][eqs1][[1]]/.Global`\[Alpha][a___][d_,e_]:>With[{k=2n,orb=Global`$OrbifoldPoint,largeN=Global`$LargeN,un=Global`$UNgroup},If[$VersionNumber<12.2,ValueQ[Global`GScoeff[orb,largeN,un][a][d,e,k]],ValueQ[Global`GScoeff[orb,largeN,un][a][d,e,k],Method->"Legacy"]]],{n,0,ord/2}]},Module[{knownorders=Count[prevord,True]},If[knownorders===1+Floor[ord/2],PrintTemporary["Gram-Schmidt coefficients at order "<>ToString[ord]<>" are already computed!"];,If[knownorders===0,print1=PrintTemporary["Computing Gram-Schmidt coefficients up to order "<>ToString[ord]];
Compute\[Alpha]temp[Quiver,ord,knownorders,eqs1];
NotebookDelete[Flatten[{print1,Table[ToExpression["tempord"<>ToString[2k]],{k,0,ord/2}]}]];
PrintTemporary["Gram-Schmidt coefficients computed at order "<>ToString[ord]];,
print1=PrintTemporary["Gram-Schmidt coefficients up to order "<>ToString[2(knownorders-1)]<>" are already computed!"<>"
Computing Gram-Schmidt coefficients up to order "<>ToString[ord]];
Compute\[Alpha]temp[Quiver,ord,knownorders,eqs1];
NotebookDelete[Flatten[{print1,Table[ToExpression["tempord"<>ToString[2k]],{k,0,ord/2}]}]];
PrintTemporary["Gram-Schmidt coefficients computed at order "<>ToString[ord]];]]]]]


\[Alpha]rules[AAA_]:=AAA/.Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][a_,k_,k_][d_,e_,f_]/;k=!=1:>(Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][a,1,1][d,e,f]/.Flatten[(Table[Global`\[Lambda][1+l]->Global`\[Lambda][k-l],{l,0,a-1}]/.Global`\[Lambda][0]:>Global`\[Lambda][a]/.Global`\[Lambda][aa_]/;Negative[aa]:>Global`\[Lambda][aa+a]/.Rule[Global`\[Lambda][aa_],Global`\[Lambda][aa_]]:>{})])/.Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][a_,k1_,k2_][d_,e_,f_]/;k1=!=1:>(Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][a,1,1+If[k2-k1<0,a+k2-k1,k2-k1]][d,e,f]/.(Table[Global`\[Lambda][l]->Global`\[Lambda][l+k1-1],{l,1,a}]/.Global`\[Lambda][aa_]/;aa>a:>Global`\[Lambda][aa-a]))


ComputeOO[Quiver_Integer/;Quiver>=1,ord_Integer/;ord>=0][{Node1_Integer/;Node1>=1,dim1_List/;AllTrue[dim1,IntegerQ]&&AllTrue[dim1,Positive]},{Node2_Integer/;Node2>=1,dim2_List/;AllTrue[dim2,IntegerQ]&&AllTrue[dim2,Positive]}]:=Which[!Global`$UNgroup&&AnyTrue[dim1,#==1&],Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" is zero in SU(N)"],!Global`$UNgroup&&AnyTrue[dim2,#==1&],Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim2,","]<>"]"],"("<>ToString[Node2]<>")"],StandardForm]<>" is zero in SU(N)"],!Global`$UNgroup&&AnyTrue[dim1,#==1&]&&AnyTrue[dim2,#==1&],Print[Style["Warning",Red],": Operators "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" and "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim2,","]<>"]"],"("<>ToString[Node2]<>")"],StandardForm]<>" are zero in SU(N)"],True,If[Node1<=Quiver &&Node2<=Quiver,If[(Plus@@dim1)===(Plus@@dim2) && Length[dim1]===Length[dim2],Module[{jj=((Plus@@XX[Quiver,Node1,dim1])(Plus@@XX[Quiver,Node2,dim2])//Expand)/.\[Chi][p_,L_]^2:>c[\[Chi][p,L],\[Chi][p,L]]/.\[Chi][p_,L1_]\[Chi][q_,L2_]:>c[\[Chi][p,L1],\[Chi][q,L2]]},ComputeCC[Quiver,ord][jj];
eqs=GSeqs[Quiver,1,dim1];
Global`old$LargeN=Global`$LargeN;
Global`$LargeN=False;
eqs2=GSeqs[Quiver,1,dim1];
Global`$LargeN=Global`old$LargeN;
Complement[eqs2//which\[Alpha],eqs//which\[Alpha]]/.Global`\[Alpha][a___][b_,c_]:>Table[(Global`GScoeff[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][a][b,c,2k]=0),{k,0,ord/2}];
Compute\[Alpha][Quiver,ord,eqs];
PrintTemporary["Summing all together..."];
fin=((jj//Ansatz[#,ord/2]&)//crules[Quiver,#]&)//\[Alpha]rules;
Global`CorrelatorOO[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup][Quiver][{Node1,dim1},{Node2,dim2}]=If[Global`$LargeN,Series[SeriesCoefficient[Normal[fin] (Global`\[Lambda][If[Global`$OrbifoldPoint,1,Node1]]^(Plus@@dim1/2) Global`\[Epsilon]^Plus@@dim1)/(8\[Pi]^2 Global`\[CapitalNu])^(Plus@@dim1/2) (Global`\[Lambda][If[Global`$OrbifoldPoint,1,Node2]]^(Plus@@dim2/2) Global`\[Epsilon]^Plus@@dim2)/(8\[Pi]^2 Global`\[CapitalNu])^(Plus@@dim2/2)//PowerExpand,{Global`\[CapitalNu],\[Infinity],0}],{Global`\[Epsilon],0,ord+Plus@@dim1+Plus@@dim2}],Series[(Collect[fin//Together,Global`\[Zeta][_],Together]//Simplify) (Global`\[Lambda][If[Global`$OrbifoldPoint,1,Node1]]^(Plus@@dim1/2) Global`\[Epsilon]^Plus@@dim1)/(8\[Pi]^2 Global`\[CapitalNu])^(Plus@@dim1/2) (Global`\[Lambda][If[Global`$OrbifoldPoint,1,Node2]]^(Plus@@dim2/2) Global`\[Epsilon]^Plus@@dim2)/(8\[Pi]^2 Global`\[CapitalNu])^(Plus@@dim2/2)//PowerExpand,{Global`\[Epsilon],0,ord+Plus@@dim1+Plus@@dim2}]]//Simplify],0],Which[Node1<= Quiver&&Node2>Quiver,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim2,","]<>"]"],"("<>ToString[Node2]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"],Node1> Quiver&&Node2<=Quiver,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"],True,Print[Style["Warning",Red],": Operators "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" and "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim2,","]<>"]"],"("<>ToString[Node2]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"]]]]


OnePFW[q_,ord_][wpos_,{a1_,a2_}]:=Module[{jj=Module[{actionexp=Series[Exp[If[EvenQ[Plus@@a2]&&Global`$LargeN,Action[q,ord/2]/.tr[X[a_]^b_.]/;OddQ[b]:>0,Action[q,ord/2]]]/.Global`\[Lambda][l_]:>Global`\[Epsilon]^2 Global`\[Lambda][l],{Global`\[Epsilon],0,ord}]//PowerExpand},If[Head[actionexp]===SeriesData,actionexp,actionexp+SeriesData[Global`\[Epsilon],0,{},0,ord+1,1]]]},((1/((jj//ExpandAll)/.tr[X[p_]^a_.]^b_.:>(\[Tau][p]@@Table[a,{i,1,b}])//.\[Tau][p_][a___]\[Tau][p_][b___]:>\[Tau][p]@@Sort[{a,b}]))((((Times@@((s[#]&/@wpos//Tally)/.{s[a_],n_}:>s[a][n][0]))(Times@@(tr[X[a1]^#]&/@a2))jj//ExpandAll)/.tr[1]->1/.tr[X[p_]^a_.]^b_.:>(\[Tau][p]@@Table[a,{i,1,b}])//.\[Tau][p_][a___]\[Tau][p_][b___]:>\[Tau][p]@@Sort[{a,b}])//.s[wlnode_][n_][0]\[Tau][wlnode_][a___]:>(s[wlnode][n]@@Sort[{a}])))/.If[Global`$UNgroup,{},\[Tau][p_][a___]/;MemberQ[List@a,1]:>0]/.\[Tau][p_][a___]/;OddQ[Total[List@a]]:>0]


\[Sigma]solve\[Lambda]exp[AA_,ord_]:=AA/.s[ii_][nn_][A___]:>Sum[(Product[1/Global`\[CapitalNu] (Global`\[Lambda][ii]/(2 Global`\[CapitalNu]))^(kk[l]/2) Global`\[Epsilon]^kk[l]/kk[l]!,{l,1,nn}]If[A===0,\[Tau][ii][Sequence@@Table[kk[l],{l,1,nn}]],\[Tau][1][A,Sequence@@Table[kk[l],{l,1,nn}]]]),Evaluate[Sequence@@Table[{kk[itemp],0,ord},{itemp,1,nn}]]]


ComputeW[Quiver_Integer/;Quiver>=1,ord_Integer/;ord>=0][posWL_List/;AllTrue[posWL,IntegerQ]&&AllTrue[posWL,Positive]]:=If[Select[posWL,#>Quiver&]==={},With[{newpos=Sort[posWL]-(Min[posWL]-1),orb=Global`$OrbifoldPoint,largeN=Global`$LargeN,un=Global`$UNgroup,texp=Global`$TranscendentalExp},Which[If[$VersionNumber<12.2,!ValueQ[Global`CorrelatorWW[orb,largeN,un,texp][Quiver][newpos]],!ValueQ[Global`CorrelatorWW[orb,largeN,un,texp][Quiver][newpos],Method->"Legacy"]],Global`CorrelatorWW[orb,largeN,un,texp][Quiver][newpos]=(If[Global`$TranscendentalExp,OnePFW[Quiver,ord][newpos,{1,{0}}]//If[Global`$UNgroup,\[Sigma]usolve,\[Sigma]solve],OnePFW[Quiver,ord][newpos,{1,{0}}]//\[Sigma]solve\[Lambda]exp[#,ord]&]//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve])/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,If[Global`$LargeN,(Coefficient[#,Global`\[CapitalNu],0]&/@c)//PowerExpand,c//PowerExpand]//Collect[#,Global`\[Zeta][__],Simplify]&,nmin,nmax,e]//Quiet,seriesOrder[Global`CorrelatorWW[orb,largeN,un,texp][Quiver][newpos]]<ord,Global`CorrelatorWW[orb,largeN,un,texp][Quiver][newpos]=Module[{WOold=Global`CorrelatorWW[orb,largeN,un,texp][Quiver][newpos]},(((If[Global`$TranscendentalExp,OnePFW[Quiver,ord][newpos,{1,{0}}]/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,Drop[Flatten[{If[Global`$LargeN,(Coefficient[#,Global`\[CapitalNu],0]&/@c),c],Table[0,{vvv,1,nmax-Length[c]}]}],seriesOrder[WOold]+1],seriesOrder[WOold]+1,nmax,e]//If[Global`$UNgroup,\[Sigma]usolve,\[Sigma]solve],(OnePFW[Quiver,ord][newpos,{1,{0}}]//\[Sigma]solve\[Lambda]exp[#,ord]&)/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,Drop[Flatten[{If[Global`$LargeN,(Coefficient[#,Global`\[CapitalNu],0]&/@c),c],Table[0,{vvv,1,nmax-Length[c]}]}],seriesOrder[WOold]+1],seriesOrder[WOold]+1,nmax,e]]//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve])/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,PowerExpand[c]//Collect[#,Global`\[Zeta][__],Simplify]&,nmin,nmax,e]//Quiet)+Normal[WOold])//Simplify],True,Series[Global`CorrelatorWW[orb,largeN,un,texp][Quiver][newpos],{Global`\[Epsilon],0,ord}]]]/.If[Global`$OrbifoldPoint,{},(Table[Global`\[Lambda][l]->Global`\[Lambda][l+Min[posWL]-1],{l,1,Quiver}]/.Global`\[Lambda][aa_]/;aa>Quiver:>Global`\[Lambda][aa-Quiver])],Print[Style["Warning",Red],": Wilson loop "<>ToString[Subscript[Global`W,"["<>StringDelete[StringDelete[ToString[Select[posWL,#>Quiver&]],"{"],"}"]<>"]"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"]]


ComputeW2[Quiver_,ord_][posWL_]:=If[Global`$LargeN,With[{newpos=Sort[posWL]-(Min[posWL]-1),orb=Global`$OrbifoldPoint,un=Global`$UNgroup,texp=Global`$TranscendentalExp},Which[If[$VersionNumber<12.2,!ValueQ[Global`CorrelatorWW2[orb,un,texp][Quiver][newpos]],!ValueQ[Global`CorrelatorWW2[orb,un,texp][Quiver][newpos],Method->"Legacy"]],Global`CorrelatorWW2[orb,un,texp][Quiver][newpos]=(If[Global`$TranscendentalExp,OnePFW[Quiver,ord][newpos,{1,{0}}]//If[Global`$UNgroup,\[Sigma]usolve,\[Sigma]solve],OnePFW[Quiver,ord][newpos,{1,{0}}]//\[Sigma]solve\[Lambda]exp[#,ord]&]//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve])/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,PowerExpand[c]//Collect[#,Global`\[Zeta][__],Simplify]&,nmin,nmax,e]//Quiet,seriesOrder[Global`CorrelatorWW2[orb,un,texp][Quiver][newpos]]<ord,Global`CorrelatorWW2[orb,un,texp][Quiver][newpos]=Module[{WOold=Global`CorrelatorWW2[orb,un,texp][Quiver][newpos]},(((If[Global`$TranscendentalExp,OnePFW[Quiver,ord][newpos,{1,{0}}]/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,Drop[Flatten[{If[Global`$LargeN,(Coefficient[#,Global`\[CapitalNu],0]&/@c),c],Table[0,{vvv,1,nmax-Length[c]}]}],seriesOrder[WOold]+1],seriesOrder[WOold]+1,nmax,e]//If[Global`$UNgroup,\[Sigma]usolve,\[Sigma]solve],(OnePFW[Quiver,ord][newpos,{1,{0}}]//\[Sigma]solve\[Lambda]exp[#,ord]&)/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,Drop[Flatten[{If[Global`$LargeN,(Coefficient[#,Global`\[CapitalNu],0]&/@c),c],Table[0,{vvv,1,nmax-Length[c]}]}],seriesOrder[WOold]+1],seriesOrder[WOold]+1,nmax,e]]//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve])/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,PowerExpand[c]//Collect[#,Global`\[Zeta][__],Simplify]&,nmin,nmax,e]//Quiet)+Normal[WOold])//Simplify],True,Series[Global`CorrelatorWW2[orb,un,texp][Quiver][newpos],{Global`\[Epsilon],0,ord}]]]/.If[Global`$OrbifoldPoint,{},(Table[Global`\[Lambda][l]->Global`\[Lambda][l+Min[posWL]-1],{l,1,Quiver}]/.Global`\[Lambda][aa_]/;aa>Quiver:>Global`\[Lambda][aa-Quiver])],ComputeW[Quiver,ord][posWL]]


NormalOrdered1PFW[Quiver_,ord_][posWL_,{posOp_,dim_}]:=If[EvenQ[Plus@@dim],(OnePFW[Quiver,ord][posWL,{posOp,dim}]-ComputeW2[Quiver,ord][posWL]OnePF[Quiver,ord][{posOp,dim}]),OnePFW[Quiver,ord][posWL,{posOp,dim}]]


ComputeWO[Quiver_Integer/;Quiver>=1,ord_Integer/;ord>=0][posWL_List/;AllTrue[posWL,IntegerQ]&&AllTrue[posWL,Positive],{posOp_Integer/;posOp>=1,dimOp_List/;AllTrue[dimOp,IntegerQ]&&AllTrue[dimOp,Positive]}]:=If[!Global`$UNgroup&&AnyTrue[dimOp,#==1&],Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dimOp,","]<>"]"],"("<>ToString[posOp]<>")"],StandardForm]<>" is zero in SU(N)"],If[Select[posWL,#>Quiver&]==={} && posOp<=Quiver,Global`old$LargeN=Global`$LargeN;
Global`$LargeN=False;
Module[{newposWL=Sort[posWL]-(Min[posWL]-1),newposOp=posOp-(Min[posWL]-1)/. 0:>Quiver/.aa_/;Negative[aa]:>Quiver+aa},Module[{jj=Plus@@XX[Quiver,newposOp,dimOp]/.\[Chi][a_,b_]:>cw[newposWL,\[Chi][a,b]]},
eqs=GSeqs[Quiver,1,dimOp];
If[EvenQ[Plus@@dimOp],printW=PrintTemporary["Computing the expectation value \[LeftAngleBracket] "<>StringReplace[ToString[Subscript["W",ToString[posWL]],StandardForm],{"{"->"[","}"->"]"}]<> " \[RightAngleBracket]"];
ComputeW2[Quiver,ord][newposWL];
NotebookDelete[printW];
PrintTemporary["Wilson loop correlators computed!"];];
ComputeCC[Quiver,ord][{eqs,jj}];
Compute\[Alpha][Quiver,ord,eqs]//Quiet;
PrintTemporary["Summing all together..."];
fin=Series[(Global`\[Lambda][If[Global`$OrbifoldPoint,1,posOp]]^(((Plus@@dimOp)-2Length[dimOp])/2)/Global`\[CapitalNu]^(((Plus@@dimOp)-2Length[dimOp])/2) Global`\[Epsilon]^((Plus@@dimOp)-2Length[dimOp])//PowerExpand)(Ansatz[jj,ord/2]//\[Alpha]rules)/.cw[pWL_,\[Chi][posO_,dim_]]:>Normal[Global`CorrSphereW[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup,Global`$TranscendentalExp][Quiver][pWL,{posO,dim}]],{Global`\[Epsilon],0,ord+(Plus@@dimOp)-2Length[dimOp]}];
Global`$LargeN=Global`old$LargeN;
Global`CorrelatorWO[Global`$OrbifoldPoint,Global`$LargeN,Global`$UNgroup,Global`$TranscendentalExp][Quiver][newposWL,{newposOp,dimOp}]=If[Global`$LargeN,Series[SeriesCoefficient[fin//Normal,{Global`\[CapitalNu],\[Infinity],0}],{Global`\[Epsilon],0,ord+(Plus@@dimOp)-2Length[dimOp]}]/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,c//Collect[#,Derivative[___][Global`W[A_]][Global`\[Lambda][B_]],Simplify]&,nmin,nmax,e]//Simplify//Quiet,fin/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,c//Collect[#,Derivative[___][Global`W[A_]][Global`\[Lambda][B_]],Simplify]&,nmin,nmax,e]//Simplify//Quiet]]]/.If[Global`$OrbifoldPoint,{},(Table[Global`\[Lambda][l]->Global`\[Lambda][l+Min[posWL]-1],{l,1,Quiver}]/.Global`\[Lambda][aa_]/;aa>Quiver:>Global`\[Lambda][aa-Quiver])],Which[Length[Select[posWL,#>Quiver&]]>0 &&  posOp<=Quiver,Print[Style["Warning",Red],": Wilson loop "<>ToString[Subscript[Global`W,"["<>StringDelete[StringDelete[ToString[Select[posWL,#>Quiver&]],"{"],"}"]<>"]"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"],Select[posWL,#>Quiver&]==={} && posOp>Quiver,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dimOp,","]<>"]"],"("<>ToString[posOp]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"],True,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dimOp,","]<>"]"],"("<>ToString[posOp]<>")"],StandardForm]<>" and Wilson loop "<>ToString[Subscript[Global`W,"["<>StringDelete[StringDelete[ToString[Select[posWL,#>Quiver&]],"{"],"}"]<>"]"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"]]]]


kernels=ParallelEvaluate[$KernelID];
SetSharedFunction[status];
SetSharedFunction[\[Tau]solve];
SetSharedFunction[\[Sigma]solve];
SetSharedFunction[\[Tau]usolve];
SetSharedFunction[\[Sigma]usolve];
DistributeDefinitions[ParallelCC];


NormalOrderedOP[Quiver_Integer/;Quiver>=1][Node1_Integer/;Node1>=1,dim1_List/;AllTrue[dim1,IntegerQ]&&AllTrue[dim1,Positive]]:=If[!Global`$UNgroup&&AnyTrue[dim1,#==1&],Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" is zero in SU(N)"],If[Node1<=Quiver,Global`old$LargeN=Global`$LargeN;
Global`$LargeN=False;
nout=(Plus@@XX[Quiver,Node1,dim1])/.\[Chi][a___]:>Global`Otilde[a];
Global`$LargeN=Global`old$LargeN;
nout,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"]]]


SphereCorrelatorOO[Quiver_Integer/;Quiver>=1,ord_Integer/;ord>=0][{Node1_Integer/;Node1>=1,dim1_List/;AllTrue[dim1,IntegerQ]&&AllTrue[dim1,Positive]},{Node2_Integer/;Node2>=1,dim2_List/;AllTrue[dim2,IntegerQ]&&AllTrue[dim2,Positive]}]:=Which[!Global`$UNgroup&&AnyTrue[dim1,#==1&],Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" is zero in SU(N)"],!Global`$UNgroup&&AnyTrue[dim2,#==1&],Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim2,","]<>"]"],"("<>ToString[Node2]<>")"],StandardForm]<>" is zero in SU(N)"],!Global`$UNgroup&&AnyTrue[dim1,#==1&]&&AnyTrue[dim2,#==1&],Print[Style["Warning",Red],": Operators "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" and "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim2,","]<>"]"],"("<>ToString[Node2]<>")"],StandardForm]<>" are zero in SU(N)"],True,If[Node1<= Quiver&&Node2<= Quiver,NormalOrdered2PF[Quiver,ord][{Node1,dim1},{Node2,dim2}]//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve],Which[Node1<= Quiver&&Node2>Quiver,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim2,","]<>"]"],"("<>ToString[Node2]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"],Node1> Quiver&&Node2<=Quiver,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"],True,Print[Style["Warning",Red],": Operators "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim1,","]<>"]"],"("<>ToString[Node1]<>")"],StandardForm]<>" and "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim2,","]<>"]"],"("<>ToString[Node2]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"]]]]


SphereCorrelatorWO[Quiver_Integer/;Quiver>=1,ord_Integer/;ord>=0][posWL_List/;AllTrue[posWL,IntegerQ]&&AllTrue[posWL,Positive],{posOP_Integer/;posOP>=1,dimOP_List/;AllTrue[dimOP,IntegerQ]&&AllTrue[dimOP,Positive]}]:=If[!Global`$UNgroup&&AnyTrue[dimOP,#==1&],Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dimOP,","]<>"]"],"("<>ToString[posOP]<>")"],StandardForm]<>" is zero in SU(N)"],If[Select[posWL,#>Quiver&]==={} && posOP<=Quiver,(If[Global`$TranscendentalExp,NormalOrdered1PFW[Quiver,ord][posWL,{posOP,dimOP}]//If[Global`$UNgroup,\[Sigma]usolve,\[Sigma]solve]//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve],(NormalOrdered1PFW[Quiver,ord][posWL,{posOP,dimOP}]//\[Sigma]solve\[Lambda]exp[#,ord]&)//If[Global`$UNgroup,\[Tau]usolve,\[Tau]solve]])/.SeriesData[a_,b_,c_,nmin_,nmax_,e_]:>SeriesData[a,b,c//If[Global`$TranscendentalExp,Collect[#,Derivative[___][Global`W[A_]][Global`\[Lambda][B_]],Together],Collect[#,Global`\[Zeta][__],Together]]&,nmin,nmax,e]//Quiet,Which[Length[Select[posWL,#>Quiver&]]>0 &&  posOP<=Quiver,Print[Style["Warning",Red],": Wilson loop "<>ToString[Subscript[Global`W,"["<>StringDelete[StringDelete[ToString[Select[posWL,#>Quiver&]],"{"],"}"]<>"]"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"],Select[posWL,#>Quiver&]==={} && posOP>Quiver,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dimOP,","]<>"]"],"("<>ToString[posOP]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"],True,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dimOP,","]<>"]"],"("<>ToString[posOP]<>")"],StandardForm]<>" and Wilson loop "<>ToString[Subscript[Global`W,"["<>StringDelete[StringDelete[ToString[Select[posWL,#>Quiver&]],"{"],"}"]<>"]"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"]]]]


GramSchmidtCoeff[Quiver_Integer/;Quiver>=1,ord_Integer/;ord>=0][Node_Integer/;Node>=1,dim_List/;AllTrue[dim,IntegerQ]&&AllTrue[dim,Positive]]:=If[!Global`$UNgroup&&AnyTrue[dim,#==1&],Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim,","]<>"]"],"("<>ToString[Node]<>")"],StandardForm]<>" is zero in SU(N)"],If[Node<=Quiver,Global`old$LargeN=Global`$LargeN;
Global`$LargeN=False;
gsout=Module[{eqs=GSeqs[Quiver,1,dim],eqs2=GSeqs[Quiver,Node,dim]},Module[{w\[Alpha]=which\[Alpha][eqs],w\[Alpha]2=which\[Alpha][eqs2]},
ComputeCC[Quiver,ord][eqs];
Compute\[Alpha][Quiver,ord,eqs];
MapThread[#1->Series[#2,{Global`\[Epsilon],0,ord}]&,{w\[Alpha]2,(w\[Alpha]2//Ansatz[#,ord/2]&)//\[Alpha]rules}]]];
Global`$LargeN=Global`old$LargeN;
gsout,Print[Style["Warning",Red],": Operator "<>ToString[Superscript[Subscript[O,"["<>StringRiffle[dim,","]<>"]"],"("<>ToString[Node]<>")"],StandardForm]<>" cannot belong to any vector multiplets in a quiver with only q="<>ToString[Quiver]<>" nodes"]]]


QUICKsaveData[file_String]:=If[!FileExistsQ[NotebookDirectory[]<>file],
Save[NotebookDirectory[]<>file,{Global`tt,Global`ss,Global`ttu,Global`ssu,Global`GScoeff,Global`CorrSphere,Global`CorrSphereW}];
Print["Saved"];,If[ChoiceDialog["File already exists. Overwrite?"],DeleteFile[NotebookDirectory[]<>file];
Save[NotebookDirectory[]<>file,{Global`tt,Global`ss,Global`ttu,Global`ssu,Global`GScoeff,Global`CorrSphere,Global`CorrSphereW}];
Print["Saved"];,$Failed]]


QUICKsaveResults[file_String]:=If[!FileExistsQ[NotebookDirectory[]<>file],
Save[NotebookDirectory[]<>file,{Global`CorrelatorOO,Global`CorrelatorWO,Global`CorrelatorWW,Global`CorrelatorWW2}];
Print["Saved"];,If[ChoiceDialog["File already exists. Overwrite?"],DeleteFile[NotebookDirectory[]<>file];
Save[NotebookDirectory[]<>file,{Global`CorrelatorOO,Global`CorrelatorWO,Global`CorrelatorWW,Global`CorrelatorWW2}];
Print["Saved"];,$Failed]]


End[]


EndPackage[]
