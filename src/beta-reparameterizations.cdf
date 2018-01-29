(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 11.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1064,         20]
NotebookDataLength[     13955,        514]
NotebookOptionsPosition[     13262,        466]
NotebookOutlinePosition[     13695,        485]
CellTagsIndexPosition[     13652,        482]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Beta", "Section"],

Cell[BoxData[
 RowBox[{
  RowBox[{"dist", "=", 
   RowBox[{"BetaDistribution", "[", 
    RowBox[{"\[Alpha]", ",", "\[Beta]"}], "]"}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[TextData[Cell[BoxData[
 FormBox[
  RowBox[{"(", 
   RowBox[{"\[Mu]", ",", 
    SuperscriptBox["\[Sigma]", "2"]}], ")"}], 
  TraditionalForm]]]], "Subsubsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{"\[Sigma]2", "==", 
       RowBox[{"Variance", "[", "dist", "]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Alpha]", ",", "\[Beta]"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{"\[Sigma]2", "==", 
       RowBox[{"Variance", "[", "dist", "]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Mu]", ",", "\[Sigma]2"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Alpha]", "\[Rule]", 
     RowBox[{"-", 
      FractionBox[
       RowBox[{"\[Mu]", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           RowBox[{"(", 
            RowBox[{
             RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}], " ", "\[Mu]"}], "+", 
          "\[Sigma]2"}], ")"}]}], "\[Sigma]2"]}]}], ",", 
    RowBox[{"\[Beta]", "\[Rule]", 
     RowBox[{
      RowBox[{"-", "1"}], "+", "\[Mu]", "+", 
      FractionBox[
       RowBox[{
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}], "2"], " ", "\[Mu]"}], 
       "\[Sigma]2"]}]}]}], "}"}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Mu]", "\[Rule]", 
     FractionBox["\[Alpha]", 
      RowBox[{"\[Alpha]", "+", "\[Beta]"}]]}], ",", 
    RowBox[{"\[Sigma]2", "\[Rule]", 
     FractionBox[
      RowBox[{"\[Alpha]", " ", "\[Beta]"}], 
      RowBox[{
       SuperscriptBox[
        RowBox[{"(", 
         RowBox[{"\[Alpha]", "+", "\[Beta]"}], ")"}], "2"], " ", 
       RowBox[{"(", 
        RowBox[{"1", "+", "\[Alpha]", "+", "\[Beta]"}], ")"}]}]]}]}], "}"}], 
  "}"}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[Cell[BoxData[
 FormBox[
  RowBox[{"(", 
   RowBox[{"\[Mu]", ",", 
    SubscriptBox["\[Mu]", "2"]}], ")"}], TraditionalForm]]]], "Subsubsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{"\[Mu]2", "\[Equal]", 
       RowBox[{"Moment", "[", 
        RowBox[{"dist", ",", "2"}], "]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Alpha]", ",", "\[Beta]"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{"\[Mu]2", "\[Equal]", 
       RowBox[{"Moment", "[", 
        RowBox[{"dist", ",", "2"}], "]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Mu]", ",", "\[Mu]2"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Alpha]", "\[Rule]", 
     FractionBox[
      RowBox[{"\[Mu]", " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "\[Mu]"}], "+", "\[Mu]2"}], ")"}]}], 
      RowBox[{
       SuperscriptBox["\[Mu]", "2"], "-", "\[Mu]2"}]]}], ",", 
    RowBox[{"\[Beta]", "\[Rule]", 
     FractionBox[
      RowBox[{
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}], " ", 
       RowBox[{"(", 
        RowBox[{"\[Mu]", "-", "\[Mu]2"}], ")"}]}], 
      RowBox[{
       SuperscriptBox["\[Mu]", "2"], "-", "\[Mu]2"}]]}]}], "}"}], 
  "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Mu]", "\[Rule]", 
     FractionBox["\[Alpha]", 
      RowBox[{"\[Alpha]", "+", "\[Beta]"}]]}], ",", 
    RowBox[{"\[Mu]2", "\[Rule]", 
     FractionBox[
      RowBox[{"\[Alpha]", " ", 
       RowBox[{"(", 
        RowBox[{"1", "+", "\[Alpha]"}], ")"}]}], 
      RowBox[{
       RowBox[{"(", 
        RowBox[{"\[Alpha]", "+", "\[Beta]"}], ")"}], " ", 
       RowBox[{"(", 
        RowBox[{"1", "+", "\[Alpha]", "+", "\[Beta]"}], ")"}]}]]}]}], "}"}], 
  "}"}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[Cell[BoxData[
 FormBox[
  RowBox[{"(", 
   RowBox[{"\[Mu]", ",", "t"}], ")"}], TraditionalForm]]]], "Subsubsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{
       RowBox[{"t", " ", 
        RowBox[{"Mean", "[", "dist", "]"}], 
        RowBox[{"(", 
         RowBox[{"1", "-", 
          RowBox[{"Mean", "[", "dist", "]"}]}], ")"}]}], "==", 
       RowBox[{"Variance", "[", "dist", "]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Alpha]", ",", "\[Beta]"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{
       RowBox[{"t", " ", 
        RowBox[{"Mean", "[", "dist", "]"}], 
        RowBox[{"(", 
         RowBox[{"1", "-", 
          RowBox[{"Mean", "[", "dist", "]"}]}], ")"}]}], "==", 
       RowBox[{"Variance", "[", "dist", "]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Mu]", ",", "t"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Alpha]", "\[Rule]", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{
        RowBox[{"-", "1"}], "+", 
        FractionBox["1", "t"]}], ")"}], " ", "\[Mu]"}]}], ",", 
    RowBox[{"\[Beta]", "\[Rule]", 
     FractionBox[
      RowBox[{
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "1"}], "+", "t"}], ")"}], " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], "t"]}]}], "}"}], 
  "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Mu]", "\[Rule]", 
     FractionBox["\[Alpha]", 
      RowBox[{"\[Alpha]", "+", "\[Beta]"}]]}], ",", 
    RowBox[{"t", "\[Rule]", 
     FractionBox["1", 
      RowBox[{"1", "+", "\[Alpha]", "+", "\[Beta]"}]]}]}], "}"}], 
  "}"}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[Cell[BoxData[
 FormBox[
  RowBox[{"(", 
   RowBox[{"\[Mu]", ",", "r"}], ")"}], TraditionalForm]]]], "Subsubsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{
       RowBox[{"r", " ", 
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{
           RowBox[{"Mean", "[", "dist", "]"}], 
           RowBox[{"(", 
            RowBox[{"1", "-", 
             RowBox[{"Mean", "[", "dist", "]"}]}], ")"}]}], ")"}], "2"]}], "==", 
       RowBox[{"Variance", "[", "dist", "]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Alpha]", ",", "\[Beta]"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{
       RowBox[{"r", " ", 
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{
           RowBox[{"Mean", "[", "dist", "]"}], 
           RowBox[{"(", 
            RowBox[{"1", "-", 
             RowBox[{"Mean", "[", "dist", "]"}]}], ")"}]}], ")"}], "2"]}], "==", 
       RowBox[{"Variance", "[", "dist", "]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Mu]", ",", "r"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Alpha]", "\[Rule]", 
     RowBox[{
      RowBox[{"-", "\[Mu]"}], "+", 
      FractionBox["1", 
       RowBox[{"r", "-", 
        RowBox[{"r", " ", "\[Mu]"}]}]]}]}], ",", 
    RowBox[{"\[Beta]", "\[Rule]", 
     RowBox[{
      RowBox[{"-", "1"}], "+", 
      FractionBox["1", 
       RowBox[{"r", " ", "\[Mu]"}]], "+", "\[Mu]"}]}]}], "}"}], 
  "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Mu]", "\[Rule]", 
     FractionBox["\[Alpha]", 
      RowBox[{"\[Alpha]", "+", "\[Beta]"}]]}], ",", 
    RowBox[{"r", "\[Rule]", 
     FractionBox[
      SuperscriptBox[
       RowBox[{"(", 
        RowBox[{"\[Alpha]", "+", "\[Beta]"}], ")"}], "2"], 
      RowBox[{"\[Alpha]", " ", "\[Beta]", " ", 
       RowBox[{"(", 
        RowBox[{"1", "+", "\[Alpha]", "+", "\[Beta]"}], ")"}]}]]}]}], "}"}], 
  "}"}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[Cell[BoxData[
 FormBox[
  RowBox[{"(", 
   RowBox[{"\[Mu]", ",", "s"}], ")"}], TraditionalForm]]]], "Subsubsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{"s", "\[Equal]", 
       RowBox[{"\[Alpha]", "+", "\[Beta]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Alpha]", ",", "\[Beta]"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{"s", "\[Equal]", 
       RowBox[{"\[Alpha]", "+", "\[Beta]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Mu]", ",", "s"}], "}"}]}], "]"}], "//", 
  "FullSimplify"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Mu]", "==", 
       RowBox[{"Mean", "[", "dist", "]"}]}], ",", 
      RowBox[{"s", "\[Equal]", 
       RowBox[{"\[Alpha]", "+", "\[Beta]"}]}], ",", 
      RowBox[{"\[Sigma]2", "\[Equal]", 
       RowBox[{"Variance", "[", "dist", "]"}]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Alpha]", ",", "\[Beta]", ",", "\[Sigma]2"}], "}"}]}], "]"}], "//",
   "FullSimplify"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Alpha]", "\[Rule]", 
     RowBox[{"s", " ", "\[Mu]"}]}], ",", 
    RowBox[{"\[Beta]", "\[Rule]", 
     RowBox[{"s", "-", 
      RowBox[{"s", " ", "\[Mu]"}]}]}]}], "}"}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Mu]", "\[Rule]", 
     FractionBox["\[Alpha]", 
      RowBox[{"\[Alpha]", "+", "\[Beta]"}]]}], ",", 
    RowBox[{"s", "\[Rule]", 
     RowBox[{"\[Alpha]", "+", "\[Beta]"}]}]}], "}"}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Alpha]", "\[Rule]", 
     RowBox[{"s", " ", "\[Mu]"}]}], ",", 
    RowBox[{"\[Beta]", "\[Rule]", 
     RowBox[{"s", "-", 
      RowBox[{"s", " ", "\[Mu]"}]}]}], ",", 
    RowBox[{"\[Sigma]2", "\[Rule]", 
     RowBox[{"-", 
      FractionBox[
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}], " ", "\[Mu]"}], 
       RowBox[{"1", "+", "s"}]]}]}]}], "}"}], "}"}]], "Output"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
},
WindowSize->{1615, 1026},
Visible->True,
ScrollingOptions->{"VerticalScrollRange"->Fit},
ShowCellBracket->Automatic,
CellContext->Notebook,
TrackCellChangeTimes->False,
FrontEndVersion->"11.0 for Linux x86 (64-bit) (September 21, 2016)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1486, 35, 23, 0, 65, "Section"],
Cell[1512, 37, 152, 4, 34, "Input"],
Cell[CellGroupData[{
Cell[1689, 45, 166, 5, 40, "Subsubsection"],
Cell[CellGroupData[{
Cell[1880, 54, 772, 24, 57, "Input"],
Cell[2655, 80, 726, 24, 54, "Output"],
Cell[3384, 106, 517, 16, 51, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[3950, 128, 158, 4, 37, "Subsubsection"],
Cell[CellGroupData[{
Cell[4133, 136, 829, 26, 57, "Input"],
Cell[4965, 164, 654, 22, 51, "Output"],
Cell[5622, 188, 541, 17, 52, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[6212, 211, 130, 3, 37, "Subsubsection"],
Cell[CellGroupData[{
Cell[6367, 218, 1096, 34, 57, "Input"],
Cell[7466, 254, 530, 19, 49, "Output"],
Cell[7999, 275, 311, 10, 51, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[8359, 291, 130, 3, 37, "Subsubsection"],
Cell[CellGroupData[{
Cell[8514, 298, 1282, 40, 71, "Input"],
Cell[9799, 340, 432, 15, 51, "Output"],
Cell[10234, 357, 489, 15, 57, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[10772, 378, 130, 3, 37, "Subsubsection"],
Cell[CellGroupData[{
Cell[10927, 385, 1251, 38, 80, "Input"],
Cell[12181, 425, 259, 8, 34, "Output"],
Cell[12443, 435, 273, 8, 48, "Output"],
Cell[12719, 445, 503, 16, 49, "Output"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature FvT0l0EHdQNfmBKblSh#7Os5 *)
