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
NotebookDataLength[    236022,       5662]
NotebookOptionsPosition[    231004,       5472]
NotebookOutlinePosition[    231602,       5497]
CellTagsIndexPosition[    231513,       5492]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Survival Distribution Simulations", "Chapter"],

Cell["\<\
This notebook numerically generates survival distributions under different \
models. These distributions are saved in the .h5 format to be imported by \
python.\
\>", "Text"],

Cell[TextData[{
 "This notebook requires the ",
 StyleBox["feature-rbsim", "Input"],
 " branch of the ",
 StyleBox["quantum-utils-mathematica", "Input"],
 " library: \n\
https://github.com/QuantumUtils/quantum-utils-mathematica/tree/feature-rbsim"
}], "Text"],

Cell[BoxData[{
 RowBox[{"Needs", "[", "\"\<QuantumUtils`\>\"", 
  "]"}], "\[IndentingNewLine]", 
 RowBox[{"Needs", "[", "\"\<RBSim`\>\"", "]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"SetDirectory", "[", 
   RowBox[{"NotebookDirectory", "[", "]"}], "]"}], ";"}]}], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"mlist", "=", 
   RowBox[{"{", 
    RowBox[{
    "1", ",", "100", ",", "200", ",", "500", ",", "1000", ",", "2000", ",", 
     "5000", ",", "10000", ",", "20000", ",", "50000"}], "}"}]}], 
  ";"}]], "Input"],

Cell[CellGroupData[{

Cell["Decay Rate Functions", "Section"],

Cell["\<\
Functions from Theorem 3 of https://arxiv.org/pdf/1703.09835.pdf to compute \
the decay rates of an RB curve for gate dependent noise.\
\>", "Text"],

Cell[BoxData[
 RowBox[{
  RowBox[{"Gu", "[", "chan_", "]"}], ":=", 
  RowBox[{"With", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"d", "=", 
      RowBox[{"InputDim", "@", "chan"}]}], "}"}], ",", 
    RowBox[{"Super", "[", 
     RowBox[{
      RowBox[{"FunctionChannel", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"chan", "[", 
          RowBox[{"#", "-", 
           RowBox[{
            RowBox[{"Tr", "[", "#", "]"}], "/", "d"}]}], "]"}], "&"}], ",", 
        RowBox[{"InputDim", "\[Rule]", "d"}]}], "]"}], ",", 
      RowBox[{"Basis", "\[Rule]", "\"\<Pauli\>\""}]}], "]"}]}], 
   "]"}]}]], "Input"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"p", "[", "gs_", "]"}], ":=", 
  RowBox[{"Abs", "@", 
   RowBox[{"First", "@", 
    RowBox[{"Eigenvalues", "[", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"First", "@", 
       RowBox[{"Super", "[", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{
          RowBox[{"Sum", "[", 
           RowBox[{
            RowBox[{
             RowBox[{"Gu", "[", 
              RowBox[{"Unitary", "@", 
               RowBox[{
                RowBox[{"GateUnitary", "[", "gs", "]"}], "[", "i", "]"}]}], 
              "]"}], "\[CircleTimes]", 
             RowBox[{"Super", "[", 
              RowBox[{
               RowBox[{
                RowBox[{"GateChannel", "[", "gs", "]"}], "[", "i", "]"}], ",", 
               RowBox[{"Basis", "\[Rule]", "\"\<Pauli\>\""}]}], "]"}]}], ",", 
            
            RowBox[{"{", 
             RowBox[{"i", ",", 
              RowBox[{"Size", "@", "gs"}]}], "}"}]}], "]"}], "/", 
          RowBox[{"Size", "[", "gs", "]"}]}], ",", "\[IndentingNewLine]", 
         RowBox[{"Basis", "\[Rule]", "\"\<Pauli\>\""}]}], 
        "\[IndentingNewLine]", "]"}]}], ",", "\[IndentingNewLine]", "1"}], 
     "\[IndentingNewLine]", "]"}]}]}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"t", "[", "gs_", "]"}], ":=", 
  RowBox[{"Abs", "@", 
   RowBox[{"First", "@", 
    RowBox[{"Eigenvalues", "[", 
     RowBox[{
      RowBox[{"First", "@", 
       RowBox[{"Super", "[", 
        RowBox[{
         RowBox[{"Sum", "[", 
          RowBox[{
           RowBox[{
            RowBox[{"GateChannel", "[", "gs", "]"}], "[", "i", "]"}], ",", 
           RowBox[{"{", 
            RowBox[{"i", ",", 
             RowBox[{"Size", "@", "gs"}]}], "}"}]}], "]"}], "/", 
         RowBox[{"Size", "[", "gs", "]"}]}], "]"}]}], ",", "1"}], 
     "]"}]}]}]}]}], "Input"],

Cell["\<\
Average gate fidelity of the depolarizing channel with p[gs] and t[gs]:\
\>", "Text"],

Cell[BoxData[
 RowBox[{
  RowBox[{"F", "[", "gs_", "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"d", "=", 
       RowBox[{"InputDim", "@", 
        RowBox[{
         RowBox[{"GateChannel", "[", "gs", "]"}], "[", "1", "]"}]}]}], ",", 
      "dep"}], "}"}], ",", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"dep", "=", 
      RowBox[{"Super", "[", 
       RowBox[{"FunctionChannel", "[", 
        RowBox[{
         RowBox[{
          RowBox[{
           RowBox[{
            RowBox[{"p", "[", "gs", "]"}], 
            RowBox[{"(", 
             RowBox[{"#", "-", 
              RowBox[{
               RowBox[{"Tr", "[", "#", "]"}], 
               RowBox[{
                RowBox[{"IdentityMatrix", "[", "d", "]"}], "/", "d"}]}]}], 
             ")"}]}], "+", " ", 
           RowBox[{
            RowBox[{"t", "[", "gs", "]"}], 
            RowBox[{"Tr", "[", "#", "]"}], 
            RowBox[{
             RowBox[{"IdentityMatrix", "[", "d", "]"}], "/", "d"}]}]}], "&"}],
          ",", 
         RowBox[{"InputDim", "\[Rule]", "d"}]}], "]"}], "]"}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"AverageGateFidelity", "[", "dep", "]"}]}]}], 
   "\[IndentingNewLine]", "]"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Solve", "[", 
  RowBox[{
   RowBox[{"p", "\[Equal]", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{
       RowBox[{"d", " ", "F"}], "-", "1"}], ")"}], "/", 
     RowBox[{"(", 
      RowBox[{"d", "-", "1"}], ")"}]}]}], ",", "F"}], "]"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{"F", "\[Rule]", 
    FractionBox[
     RowBox[{"1", "-", "p", "+", 
      RowBox[{"d", " ", "p"}]}], "d"]}], "}"}], "}"}]], "Output"]
}, Open  ]]
}, Closed]],

Cell[CellGroupData[{

Cell["Depolarizing Error Model", "Section"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"depolarizingGate", "[", "p_", "]"}], ":=", 
  RowBox[{"Super", "[", 
   RowBox[{"FunctionChannel", "[", 
    RowBox[{
     RowBox[{
      RowBox[{
       RowBox[{
        RowBox[{"(", 
         RowBox[{"1", "-", "p"}], ")"}], "#"}], "+", 
       RowBox[{"p", " ", 
        RowBox[{
         RowBox[{"Tr", "[", "#", "]"}], "/", "2"}]}]}], "&"}], ",", 
     RowBox[{"InputDim", "\[Rule]", "2"}]}], "]"}], 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateNoise", "=", 
   RowBox[{"IndependentNoise", "[", 
    RowBox[{"depolarizingGate", "[", "0.0002", "]"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"gateChannels", "=", 
    RowBox[{
     RowBox[{
      RowBox[{
       RowBox[{"Unitary", "[", 
        RowBox[{
         RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", "#",
          "]"}], "]"}], ".", 
       RowBox[{"First", "[", "gateNoise", "]"}]}], "&"}], "/@", 
     RowBox[{"Range", "[", "12", "]"}]}]}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gs", "=", 
   RowBox[{"GateNoiseGateSet", "[", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"GateSetName", "\[Rule]", "\"\<Depolarizing Noise\>\""}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"Size", "\[Rule]", "12"}], ",", "\[IndentingNewLine]", 
     RowBox[{"Dimension", "\[Rule]", "2"}], ",", "\[IndentingNewLine]", 
     RowBox[{"GateProduct", "\[Rule]", 
      RowBox[{"GateProduct", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateInverse", "\[Rule]", 
      RowBox[{"GateInverse", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateUnitary", "\[Rule]", 
      RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateNoise", "\[Rule]", "gateNoise"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateChannel", "\[Rule]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"gateChannels", "[", 
         RowBox[{"[", 
          RowBox[{"Last", "[", 
           RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}], ",",
      "\[IndentingNewLine]", 
     RowBox[{"GateNoiseMemoryDepth", "\[Rule]", "1"}]}], 
    "\[IndentingNewLine]", "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"{", 
  RowBox[{
   RowBox[{"p", "[", "gs", "]"}], ",", 
   RowBox[{"t", "[", "gs", "]"}], ",", 
   RowBox[{"1", "-", 
    RowBox[{"F", "[", "gs", "]"}]}]}], "}"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0.9997999999999996`", ",", "0.9999999999999999`", ",", 
   "0.00010000000000021103`"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"protocol", "=", 
   RowBox[{"RBProtocol", "[", 
    RowBox[{"mlist", ",", "5000", ",", "1", ",", 
     RowBox[{"TP", "[", "U", "]"}], ",", 
     RowBox[{"0.99", 
      RowBox[{"TP", "[", "U", "]"}]}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"survivals", "=", 
   RowBox[{"SimulateProtocol", "[", 
    RowBox[{"gs", ",", "protocol", ",", 
     RowBox[{
     "SimulationExportName", "\[Rule]", 
      "\"\<../data/depolarizing_model_survivals\>\""}]}], "]"}]}], 
  ";"}]}], "Input"],

Cell[BoxData["\<\"Saved ../data/depolarizing_model_survivals.h5\"\>"], "Print"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"Histogram", "[", 
    RowBox[{"#", ",", 
     RowBox[{"{", 
      RowBox[{"0", ",", "1", ",", "0.01"}], "}"}], ",", 
     RowBox[{"ImageSize", "\[Rule]", "400"}]}], "]"}], "&"}], "/@", 
  RowBox[{"survivals", "[", 
   RowBox[{"[", 
    RowBox[{"All", ",", "1", ",", "All", ",", "1"}], "]"}], "]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.98, 0}, {0.99, 5000},
        RoundingRadius->0]}, {}, {}}, {}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.97, 0}, {0.98, 2029},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 2971},
        RoundingRadius->0]}, {}, {}}, {{}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.96, 0}, {0.97, 1300},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 3700},
        RoundingRadius->0]}, {}, {}}, {{}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.93, 0}, {0.9400000000000001, 51},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 4949},
        RoundingRadius->0]}, {}, {}}, {{}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.89, 0}, {0.9, 2204},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 2796},
        RoundingRadius->0]}, {}, {}}, {{}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.81, 0}, {0.8200000000000001, 3},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 4754},
        RoundingRadius->0], 
       RectangleBox[{0.8300000000000001, 0}, {0.84, 243},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.66, 0}, {0.67, 18},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 4322},
        RoundingRadius->0], 
       RectangleBox[{0.68, 0}, {0.6900000000000001, 660},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.55, 0}, {0.56, 1184},
        RoundingRadius->0], 
       RectangleBox[{0.56, 0}, {0.5700000000000001, 3804},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 12},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.49, 0}, {0.5, 381},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 4532},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 87},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.48, 0}, {0.49, 202},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 4592},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 206},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}]}], "}"}]], "Output"]
}, Open  ]]
}, Closed]],

Cell[CellGroupData[{

Cell["Overrotation Error Model", "Section"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"gateChannels", "=", 
   RowBox[{"Table", "[", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"If", "[", 
      RowBox[{
       RowBox[{"DiagonalMatrixQ", "[", 
        RowBox[{
         RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", "i",
          "]"}], "]"}], ",", "\[IndentingNewLine]", 
       RowBox[{"Super", "@", 
        RowBox[{"Unitary", "@", 
         RowBox[{
          RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", 
          "i", "]"}]}]}], ",", "\[IndentingNewLine]", 
       RowBox[{"Super", "@", 
        RowBox[{"Unitary", "[", 
         RowBox[{"MatrixPower", "[", 
          RowBox[{
           RowBox[{
            RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", 
            "i", "]"}], ",", "1.011132"}], "]"}], "]"}]}]}], 
      "\[IndentingNewLine]", "]"}], ",", "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"i", ",", "12"}], "}"}]}], "\[IndentingNewLine]", "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateNoise", "=", 
   RowBox[{"Table", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"Unitary", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", "i",
          "]"}], "\[ConjugateTranspose]"}], "]"}], ".", 
      RowBox[{"gateChannels", "[", 
       RowBox[{"[", "i", "]"}], "]"}]}], ",", 
     RowBox[{"{", 
      RowBox[{"i", ",", "12"}], "}"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gs", "=", 
   RowBox[{"GateNoiseGateSet", "[", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"GateSetName", "\[Rule]", "\"\<Overrotation Noise\>\""}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"Size", "\[Rule]", "12"}], ",", "\[IndentingNewLine]", 
     RowBox[{"Dimension", "\[Rule]", "2"}], ",", "\[IndentingNewLine]", 
     RowBox[{"GateProduct", "\[Rule]", 
      RowBox[{"GateProduct", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateInverse", "\[Rule]", 
      RowBox[{"GateInverse", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateUnitary", "\[Rule]", 
      RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateNoise", "\[Rule]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"gateNoise", "[", 
         RowBox[{"[", 
          RowBox[{"Last", "[", 
           RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}], ",",
      "\[IndentingNewLine]", 
     RowBox[{"GateChannel", "\[Rule]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"gateChannels", "[", 
         RowBox[{"[", 
          RowBox[{"Last", "[", 
           RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}], ",",
      "\[IndentingNewLine]", 
     RowBox[{"GateNoiseMemoryDepth", "\[Rule]", "1"}]}], 
    "\[IndentingNewLine]", "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"{", 
  RowBox[{
   RowBox[{"p", "[", "gs", "]"}], ",", 
   RowBox[{"t", "[", "gs", "]"}], ",", 
   RowBox[{"1", "-", 
    RowBox[{"F", "[", "gs", "]"}]}]}], "}"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0.9997999990767658`", ",", "0.9999999999999999`", ",", 
   "0.0001000004616171779`"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"protocol", "=", 
   RowBox[{"RBProtocol", "[", 
    RowBox[{"mlist", ",", "5000", ",", "1", ",", 
     RowBox[{"TP", "[", "U", "]"}], ",", 
     RowBox[{"0.99", 
      RowBox[{"TP", "[", "U", "]"}]}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"survivals", "=", 
   RowBox[{"SimulateProtocol", "[", 
    RowBox[{"gs", ",", "protocol", ",", 
     RowBox[{
     "SimulationExportName", "\[Rule]", 
      "\"\<../data/overrotation_model_survivals\>\""}]}], "]"}]}], 
  ";"}]}], "Input"],

Cell[BoxData["\<\"Saved ../data/overrotation_model_survivals.h5\"\>"], "Print"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"Histogram", "[", 
    RowBox[{"#", ",", 
     RowBox[{"{", 
      RowBox[{"0", ",", "1", ",", "0.01"}], "}"}], ",", 
     RowBox[{"ImageSize", "\[Rule]", "400"}]}], "]"}], "&"}], "/@", 
  RowBox[{"survivals", "[", 
   RowBox[{"[", 
    RowBox[{"All", ",", "1", ",", "All", ",", "1"}], "]"}], "]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.98, 0}, {0.99, 820},
        RoundingRadius->0], RectangleBox[{0.99, 0}, {1., 4180},
        RoundingRadius->0]}, {}, {}}, {{}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.9, 0}, {0.91, 1},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 1},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 8},
        RoundingRadius->0], RectangleBox[{0.93, 0}, {0.9400000000000001, 23},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 67},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 132},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 428},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 1153},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 3187},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.81, 0}, {0.8200000000000001, 1},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 3},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 1},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 3},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 2},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 4},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 8},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 11},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 17},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 30},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 65},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 94},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 168},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 271},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 483},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 741},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 1195},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 1903},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.5700000000000001, 0}, {0.58, 1},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 1},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 1},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 2},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 1},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 3},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 4},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 2},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 6},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 6},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 6},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 11},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 5},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 14},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 5},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 16},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 23},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 24},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 32},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 48},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 45},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 63},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 72},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 81},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 130},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 156},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 191},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 214},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 295},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 336},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 383},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 520},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 616},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 741},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 946},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.36, 0}, {0.37, 1},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 1},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 1},
        RoundingRadius->0], RectangleBox[{0.47000000000000003, 0}, {0.48, 1},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 1},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 2},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 3},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 3},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 5},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 4},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 3},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 5},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 4},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 6},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 1},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 4},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 6},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 10},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 14},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 15},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 13},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 11},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 18},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 19},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 22},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 27},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 27},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 28},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 42},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 40},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 47},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 48},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 45},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 66},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 58},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 70},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 77},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 79},
        RoundingRadius->0], 
       RectangleBox[{0.8300000000000001, 0}, {0.84, 133},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 111},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 117},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 136},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 161},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 181},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 219},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 243},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 251},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 264},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 320},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 323},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 370},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 416},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 450},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 478},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.11, 0}, {0.12, 1},
        RoundingRadius->0], RectangleBox[{0.13, 0}, {0.14, 1},
        RoundingRadius->0], RectangleBox[{0.15, 0}, {0.16, 1},
        RoundingRadius->0], RectangleBox[{0.16, 0}, {0.17, 1},
        RoundingRadius->0], RectangleBox[{0.18, 0}, {0.19, 1},
        RoundingRadius->0], RectangleBox[{0.19, 0}, {0.2, 1},
        RoundingRadius->0], RectangleBox[{0.2, 0}, {0.21, 1},
        RoundingRadius->0], RectangleBox[{0.22, 0}, {0.23, 4},
        RoundingRadius->0], RectangleBox[{0.23, 0}, {0.24, 2},
        RoundingRadius->0], RectangleBox[{0.24, 0}, {0.25, 3},
        RoundingRadius->0], RectangleBox[{0.25, 0}, {0.26, 3},
        RoundingRadius->0], RectangleBox[{0.26, 0}, {0.27, 4},
        RoundingRadius->0], RectangleBox[{0.27, 0}, {0.28, 5},
        RoundingRadius->0], RectangleBox[{0.28, 0}, {0.29, 2},
        RoundingRadius->0], RectangleBox[{0.29, 0}, {0.3, 5},
        RoundingRadius->0], RectangleBox[{0.3, 0}, {0.31, 1},
        RoundingRadius->0], RectangleBox[{0.31, 0}, {0.32, 2},
        RoundingRadius->0], RectangleBox[{0.32, 0}, {0.33, 2},
        RoundingRadius->0], RectangleBox[{0.33, 0}, {0.34, 5},
        RoundingRadius->0], RectangleBox[{0.34, 0}, {0.35000000000000003, 4},
        RoundingRadius->0], RectangleBox[{0.35000000000000003, 0}, {0.36, 3},
        RoundingRadius->0], RectangleBox[{0.36, 0}, {0.37, 10},
        RoundingRadius->0], RectangleBox[{0.37, 0}, {0.38, 2},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 9},
        RoundingRadius->0], RectangleBox[{0.39, 0}, {0.4, 12},
        RoundingRadius->0], RectangleBox[{0.4, 0}, {0.41000000000000003, 9},
        RoundingRadius->0], RectangleBox[{0.41000000000000003, 0}, {0.42, 8},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 9},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 5},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 8},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 10},
        RoundingRadius->0], 
       RectangleBox[{0.46, 0}, {0.47000000000000003, 13},
        RoundingRadius->0], 
       RectangleBox[{0.47000000000000003, 0}, {0.48, 13},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 14},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 12},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 13},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 15},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 18},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 23},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 18},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 21},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 15},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 19},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 26},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 28},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 30},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 36},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 39},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 53},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 41},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 44},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 41},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 54},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 64},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 60},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 64},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 47},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 61},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 89},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 77},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 79},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 100},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 84},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 102},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 99},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 100},
        RoundingRadius->0], 
       RectangleBox[{0.81, 0}, {0.8200000000000001, 130},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 120},
        RoundingRadius->0], 
       RectangleBox[{0.8300000000000001, 0}, {0.84, 130},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 122},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 135},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 129},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 142},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 175},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 155},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 188},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 201},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 186},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 217},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 240},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 249},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 237},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 243},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 255},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0., 0}, {0.01, 2},
        RoundingRadius->0], RectangleBox[{0.01, 0}, {0.02, 11},
        RoundingRadius->0], RectangleBox[{0.02, 0}, {0.03, 5},
        RoundingRadius->0], RectangleBox[{0.03, 0}, {0.04, 5},
        RoundingRadius->0], RectangleBox[{0.04, 0}, {0.05, 9},
        RoundingRadius->0], RectangleBox[{0.05, 0}, {0.06, 9},
        RoundingRadius->0], RectangleBox[{0.06, 0}, {0.07, 9},
        RoundingRadius->0], RectangleBox[{0.07, 0}, {0.08, 14},
        RoundingRadius->0], RectangleBox[{0.08, 0}, {0.09, 9},
        RoundingRadius->0], RectangleBox[{0.09, 0}, {0.1, 8},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 9},
        RoundingRadius->0], RectangleBox[{0.11, 0}, {0.12, 10},
        RoundingRadius->0], RectangleBox[{0.12, 0}, {0.13, 9},
        RoundingRadius->0], RectangleBox[{0.13, 0}, {0.14, 21},
        RoundingRadius->0], RectangleBox[{0.14, 0}, {0.15, 13},
        RoundingRadius->0], RectangleBox[{0.15, 0}, {0.16, 13},
        RoundingRadius->0], RectangleBox[{0.16, 0}, {0.17, 12},
        RoundingRadius->0], RectangleBox[{0.17, 0}, {0.18, 17},
        RoundingRadius->0], RectangleBox[{0.18, 0}, {0.19, 11},
        RoundingRadius->0], RectangleBox[{0.19, 0}, {0.2, 14},
        RoundingRadius->0], RectangleBox[{0.2, 0}, {0.21, 26},
        RoundingRadius->0], RectangleBox[{0.21, 0}, {0.22, 10},
        RoundingRadius->0], RectangleBox[{0.22, 0}, {0.23, 26},
        RoundingRadius->0], RectangleBox[{0.23, 0}, {0.24, 18},
        RoundingRadius->0], RectangleBox[{0.24, 0}, {0.25, 24},
        RoundingRadius->0], RectangleBox[{0.25, 0}, {0.26, 25},
        RoundingRadius->0], RectangleBox[{0.26, 0}, {0.27, 19},
        RoundingRadius->0], RectangleBox[{0.27, 0}, {0.28, 26},
        RoundingRadius->0], RectangleBox[{0.28, 0}, {0.29, 20},
        RoundingRadius->0], RectangleBox[{0.29, 0}, {0.3, 23},
        RoundingRadius->0], RectangleBox[{0.3, 0}, {0.31, 25},
        RoundingRadius->0], RectangleBox[{0.31, 0}, {0.32, 33},
        RoundingRadius->0], RectangleBox[{0.32, 0}, {0.33, 30},
        RoundingRadius->0], RectangleBox[{0.33, 0}, {0.34, 30},
        RoundingRadius->0], 
       RectangleBox[{0.34, 0}, {0.35000000000000003, 27},
        RoundingRadius->0], 
       RectangleBox[{0.35000000000000003, 0}, {0.36, 33},
        RoundingRadius->0], RectangleBox[{0.36, 0}, {0.37, 38},
        RoundingRadius->0], RectangleBox[{0.37, 0}, {0.38, 39},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 39},
        RoundingRadius->0], RectangleBox[{0.39, 0}, {0.4, 46},
        RoundingRadius->0], RectangleBox[{0.4, 0}, {0.41000000000000003, 23},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 31},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 31},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 38},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 36},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 29},
        RoundingRadius->0], 
       RectangleBox[{0.46, 0}, {0.47000000000000003, 54},
        RoundingRadius->0], 
       RectangleBox[{0.47000000000000003, 0}, {0.48, 43},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 58},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 38},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 55},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 42},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 43},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 41},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 44},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 66},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 59},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 43},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 50},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 75},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 62},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 51},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 57},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 62},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 68},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 58},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 53},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 77},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 76},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 67},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 67},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 74},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 80},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 75},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 67},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 65},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 80},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 77},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 74},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 78},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 90},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 92},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 99},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 86},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 106},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 76},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 112},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 101},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 91},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 104},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 101},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 98},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 117},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 111},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 114},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 112},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 124},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 109},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 123},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0., 0}, {0.01, 29},
        RoundingRadius->0], RectangleBox[{0.01, 0}, {0.02, 23},
        RoundingRadius->0], RectangleBox[{0.02, 0}, {0.03, 24},
        RoundingRadius->0], RectangleBox[{0.03, 0}, {0.04, 39},
        RoundingRadius->0], RectangleBox[{0.04, 0}, {0.05, 35},
        RoundingRadius->0], RectangleBox[{0.05, 0}, {0.06, 32},
        RoundingRadius->0], RectangleBox[{0.06, 0}, {0.07, 33},
        RoundingRadius->0], RectangleBox[{0.07, 0}, {0.08, 31},
        RoundingRadius->0], RectangleBox[{0.08, 0}, {0.09, 34},
        RoundingRadius->0], RectangleBox[{0.09, 0}, {0.1, 29},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 45},
        RoundingRadius->0], RectangleBox[{0.11, 0}, {0.12, 32},
        RoundingRadius->0], RectangleBox[{0.12, 0}, {0.13, 38},
        RoundingRadius->0], RectangleBox[{0.13, 0}, {0.14, 34},
        RoundingRadius->0], RectangleBox[{0.14, 0}, {0.15, 45},
        RoundingRadius->0], RectangleBox[{0.15, 0}, {0.16, 40},
        RoundingRadius->0], RectangleBox[{0.16, 0}, {0.17, 40},
        RoundingRadius->0], RectangleBox[{0.17, 0}, {0.18, 41},
        RoundingRadius->0], RectangleBox[{0.18, 0}, {0.19, 42},
        RoundingRadius->0], RectangleBox[{0.19, 0}, {0.2, 33},
        RoundingRadius->0], RectangleBox[{0.2, 0}, {0.21, 27},
        RoundingRadius->0], RectangleBox[{0.21, 0}, {0.22, 43},
        RoundingRadius->0], RectangleBox[{0.22, 0}, {0.23, 51},
        RoundingRadius->0], RectangleBox[{0.23, 0}, {0.24, 47},
        RoundingRadius->0], RectangleBox[{0.24, 0}, {0.25, 41},
        RoundingRadius->0], RectangleBox[{0.25, 0}, {0.26, 38},
        RoundingRadius->0], RectangleBox[{0.26, 0}, {0.27, 32},
        RoundingRadius->0], RectangleBox[{0.27, 0}, {0.28, 34},
        RoundingRadius->0], RectangleBox[{0.28, 0}, {0.29, 43},
        RoundingRadius->0], RectangleBox[{0.29, 0}, {0.3, 39},
        RoundingRadius->0], RectangleBox[{0.3, 0}, {0.31, 43},
        RoundingRadius->0], RectangleBox[{0.31, 0}, {0.32, 44},
        RoundingRadius->0], RectangleBox[{0.32, 0}, {0.33, 41},
        RoundingRadius->0], RectangleBox[{0.33, 0}, {0.34, 37},
        RoundingRadius->0], 
       RectangleBox[{0.34, 0}, {0.35000000000000003, 40},
        RoundingRadius->0], 
       RectangleBox[{0.35000000000000003, 0}, {0.36, 46},
        RoundingRadius->0], RectangleBox[{0.36, 0}, {0.37, 45},
        RoundingRadius->0], RectangleBox[{0.37, 0}, {0.38, 42},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 52},
        RoundingRadius->0], RectangleBox[{0.39, 0}, {0.4, 44},
        RoundingRadius->0], RectangleBox[{0.4, 0}, {0.41000000000000003, 45},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 57},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 48},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 49},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 59},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 40},
        RoundingRadius->0], 
       RectangleBox[{0.46, 0}, {0.47000000000000003, 60},
        RoundingRadius->0], 
       RectangleBox[{0.47000000000000003, 0}, {0.48, 51},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 46},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 53},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 53},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 46},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 53},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 47},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 42},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 56},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 40},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 61},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 53},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 60},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 61},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 46},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 58},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 60},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 51},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 62},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 64},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 65},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 66},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 66},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 60},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 49},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 64},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 60},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 67},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 59},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 47},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 58},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 64},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 64},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 75},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 75},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 65},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 61},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 65},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 61},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 63},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 61},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 58},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 69},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 68},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 75},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 72},
        RoundingRadius->0], RectangleBox[{0.93, 0}, {0.9400000000000001, 61},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 71},
        RoundingRadius->0], RectangleBox[{0.9500000000000001, 0}, {0.96, 53},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 67},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 77},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 65},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0., 0}, {0.01, 50},
        RoundingRadius->0], RectangleBox[{0.01, 0}, {0.02, 39},
        RoundingRadius->0], RectangleBox[{0.02, 0}, {0.03, 56},
        RoundingRadius->0], RectangleBox[{0.03, 0}, {0.04, 41},
        RoundingRadius->0], RectangleBox[{0.04, 0}, {0.05, 51},
        RoundingRadius->0], RectangleBox[{0.05, 0}, {0.06, 44},
        RoundingRadius->0], RectangleBox[{0.06, 0}, {0.07, 49},
        RoundingRadius->0], RectangleBox[{0.07, 0}, {0.08, 39},
        RoundingRadius->0], RectangleBox[{0.08, 0}, {0.09, 44},
        RoundingRadius->0], RectangleBox[{0.09, 0}, {0.1, 44},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 49},
        RoundingRadius->0], RectangleBox[{0.11, 0}, {0.12, 45},
        RoundingRadius->0], RectangleBox[{0.12, 0}, {0.13, 61},
        RoundingRadius->0], RectangleBox[{0.13, 0}, {0.14, 35},
        RoundingRadius->0], RectangleBox[{0.14, 0}, {0.15, 52},
        RoundingRadius->0], RectangleBox[{0.15, 0}, {0.16, 60},
        RoundingRadius->0], RectangleBox[{0.16, 0}, {0.17, 48},
        RoundingRadius->0], RectangleBox[{0.17, 0}, {0.18, 46},
        RoundingRadius->0], RectangleBox[{0.18, 0}, {0.19, 56},
        RoundingRadius->0], RectangleBox[{0.19, 0}, {0.2, 52},
        RoundingRadius->0], RectangleBox[{0.2, 0}, {0.21, 48},
        RoundingRadius->0], RectangleBox[{0.21, 0}, {0.22, 48},
        RoundingRadius->0], RectangleBox[{0.22, 0}, {0.23, 58},
        RoundingRadius->0], RectangleBox[{0.23, 0}, {0.24, 56},
        RoundingRadius->0], RectangleBox[{0.24, 0}, {0.25, 51},
        RoundingRadius->0], RectangleBox[{0.25, 0}, {0.26, 55},
        RoundingRadius->0], RectangleBox[{0.26, 0}, {0.27, 54},
        RoundingRadius->0], RectangleBox[{0.27, 0}, {0.28, 46},
        RoundingRadius->0], RectangleBox[{0.28, 0}, {0.29, 42},
        RoundingRadius->0], RectangleBox[{0.29, 0}, {0.3, 49},
        RoundingRadius->0], RectangleBox[{0.3, 0}, {0.31, 61},
        RoundingRadius->0], RectangleBox[{0.31, 0}, {0.32, 52},
        RoundingRadius->0], RectangleBox[{0.32, 0}, {0.33, 49},
        RoundingRadius->0], RectangleBox[{0.33, 0}, {0.34, 48},
        RoundingRadius->0], 
       RectangleBox[{0.34, 0}, {0.35000000000000003, 42},
        RoundingRadius->0], 
       RectangleBox[{0.35000000000000003, 0}, {0.36, 55},
        RoundingRadius->0], RectangleBox[{0.36, 0}, {0.37, 50},
        RoundingRadius->0], RectangleBox[{0.37, 0}, {0.38, 59},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 53},
        RoundingRadius->0], RectangleBox[{0.39, 0}, {0.4, 53},
        RoundingRadius->0], RectangleBox[{0.4, 0}, {0.41000000000000003, 48},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 42},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 43},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 46},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 56},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 53},
        RoundingRadius->0], 
       RectangleBox[{0.46, 0}, {0.47000000000000003, 53},
        RoundingRadius->0], 
       RectangleBox[{0.47000000000000003, 0}, {0.48, 55},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 58},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 43},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 50},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 51},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 42},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 60},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 46},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 51},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 55},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 47},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 61},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 49},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 54},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 53},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 59},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 57},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 43},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 52},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 46},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 43},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 55},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 52},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 58},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 53},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 55},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 45},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 49},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 39},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 52},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 49},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 46},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 49},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 54},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 50},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 51},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 66},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 45},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 55},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 50},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 51},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 46},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 58},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 66},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 56},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 63},
        RoundingRadius->0], RectangleBox[{0.93, 0}, {0.9400000000000001, 44},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 58},
        RoundingRadius->0], RectangleBox[{0.9500000000000001, 0}, {0.96, 57},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 40},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 42},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 40},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0., 0}, {0.01, 46},
        RoundingRadius->0], RectangleBox[{0.01, 0}, {0.02, 45},
        RoundingRadius->0], RectangleBox[{0.02, 0}, {0.03, 45},
        RoundingRadius->0], RectangleBox[{0.03, 0}, {0.04, 46},
        RoundingRadius->0], RectangleBox[{0.04, 0}, {0.05, 46},
        RoundingRadius->0], RectangleBox[{0.05, 0}, {0.06, 58},
        RoundingRadius->0], RectangleBox[{0.06, 0}, {0.07, 50},
        RoundingRadius->0], RectangleBox[{0.07, 0}, {0.08, 57},
        RoundingRadius->0], RectangleBox[{0.08, 0}, {0.09, 60},
        RoundingRadius->0], RectangleBox[{0.09, 0}, {0.1, 50},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 46},
        RoundingRadius->0], RectangleBox[{0.11, 0}, {0.12, 47},
        RoundingRadius->0], RectangleBox[{0.12, 0}, {0.13, 46},
        RoundingRadius->0], RectangleBox[{0.13, 0}, {0.14, 61},
        RoundingRadius->0], RectangleBox[{0.14, 0}, {0.15, 59},
        RoundingRadius->0], RectangleBox[{0.15, 0}, {0.16, 54},
        RoundingRadius->0], RectangleBox[{0.16, 0}, {0.17, 44},
        RoundingRadius->0], RectangleBox[{0.17, 0}, {0.18, 55},
        RoundingRadius->0], RectangleBox[{0.18, 0}, {0.19, 49},
        RoundingRadius->0], RectangleBox[{0.19, 0}, {0.2, 49},
        RoundingRadius->0], RectangleBox[{0.2, 0}, {0.21, 50},
        RoundingRadius->0], RectangleBox[{0.21, 0}, {0.22, 61},
        RoundingRadius->0], RectangleBox[{0.22, 0}, {0.23, 40},
        RoundingRadius->0], RectangleBox[{0.23, 0}, {0.24, 40},
        RoundingRadius->0], RectangleBox[{0.24, 0}, {0.25, 43},
        RoundingRadius->0], RectangleBox[{0.25, 0}, {0.26, 55},
        RoundingRadius->0], RectangleBox[{0.26, 0}, {0.27, 57},
        RoundingRadius->0], RectangleBox[{0.27, 0}, {0.28, 56},
        RoundingRadius->0], RectangleBox[{0.28, 0}, {0.29, 45},
        RoundingRadius->0], RectangleBox[{0.29, 0}, {0.3, 44},
        RoundingRadius->0], RectangleBox[{0.3, 0}, {0.31, 43},
        RoundingRadius->0], RectangleBox[{0.31, 0}, {0.32, 63},
        RoundingRadius->0], RectangleBox[{0.32, 0}, {0.33, 63},
        RoundingRadius->0], RectangleBox[{0.33, 0}, {0.34, 39},
        RoundingRadius->0], 
       RectangleBox[{0.34, 0}, {0.35000000000000003, 54},
        RoundingRadius->0], 
       RectangleBox[{0.35000000000000003, 0}, {0.36, 51},
        RoundingRadius->0], RectangleBox[{0.36, 0}, {0.37, 62},
        RoundingRadius->0], RectangleBox[{0.37, 0}, {0.38, 50},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 49},
        RoundingRadius->0], RectangleBox[{0.39, 0}, {0.4, 54},
        RoundingRadius->0], RectangleBox[{0.4, 0}, {0.41000000000000003, 52},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 54},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 56},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 48},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 55},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 46},
        RoundingRadius->0], 
       RectangleBox[{0.46, 0}, {0.47000000000000003, 50},
        RoundingRadius->0], 
       RectangleBox[{0.47000000000000003, 0}, {0.48, 51},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 63},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 58},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 59},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 38},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 45},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 54},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 53},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 46},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 54},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 47},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 52},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 53},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 48},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 38},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 63},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 52},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 41},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 43},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 50},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 54},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 46},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 45},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 43},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 48},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 56},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 53},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 44},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 47},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 50},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 61},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 50},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 44},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 47},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 43},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 57},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 47},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 68},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 47},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 52},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 51},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 50},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 43},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 49},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 57},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 58},
        RoundingRadius->0], RectangleBox[{0.93, 0}, {0.9400000000000001, 51},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 46},
        RoundingRadius->0], RectangleBox[{0.9500000000000001, 0}, {0.96, 45},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 48},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 54},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 45},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}]}], "}"}]], "Output"]
}, Open  ]]
}, Closed]],

Cell[CellGroupData[{

Cell["Overrotation with Dephasing Error Model", "Section"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"dephasingGate", "[", "p_", "]"}], ":=", 
  RowBox[{"Kraus", "[", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{
      RowBox[{"Sqrt", "[", 
       RowBox[{"1", "-", "p"}], "]"}], 
      RowBox[{"TP", "[", "I", "]"}]}], ",", 
     RowBox[{
      RowBox[{"Sqrt", "[", "p", "]"}], 
      RowBox[{"TP", "[", "Z", "]"}]}]}], "}"}], 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateChannels", "=", 
   RowBox[{"Table", "[", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
      RowBox[{"If", "[", 
       RowBox[{
        RowBox[{"DiagonalMatrixQ", "[", 
         RowBox[{
          RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", 
          "i", "]"}], "]"}], ",", "\[IndentingNewLine]", 
        RowBox[{"Super", "@", 
         RowBox[{"Unitary", "@", 
          RowBox[{
           RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", 
           "i", "]"}]}]}], ",", "\[IndentingNewLine]", 
        RowBox[{"Super", "@", 
         RowBox[{"Unitary", "[", 
          RowBox[{"MatrixPower", "[", 
           RowBox[{
            RowBox[{
             RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", 
             "i", "]"}], ",", "1.01"}], "]"}], "]"}]}]}], 
       "\[IndentingNewLine]", "]"}], ".", 
      RowBox[{"dephasingGate", "[", "0.000028954", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"i", ",", "12"}], "}"}]}], "\[IndentingNewLine]", "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"gateNoise", "=", 
    RowBox[{"Table", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"Unitary", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", 
          "i", "]"}], "\[ConjugateTranspose]"}], "]"}], ".", 
       RowBox[{"gateChannels", "[", 
        RowBox[{"[", "i", "]"}], "]"}]}], ",", 
      RowBox[{"{", 
       RowBox[{"i", ",", "12"}], "}"}]}], "]"}]}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gs", "=", 
   RowBox[{"GateNoiseGateSet", "[", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
     "GateSetName", "\[Rule]", "\"\<Overrotation+Dephasing Noise\>\""}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"Size", "\[Rule]", "12"}], ",", "\[IndentingNewLine]", 
     RowBox[{"Dimension", "\[Rule]", "2"}], ",", "\[IndentingNewLine]", 
     RowBox[{"GateProduct", "\[Rule]", 
      RowBox[{"GateProduct", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateInverse", "\[Rule]", 
      RowBox[{"GateInverse", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateUnitary", "\[Rule]", 
      RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateNoise", "\[Rule]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"gateNoise", "[", 
         RowBox[{"[", 
          RowBox[{"Last", "[", 
           RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}], ",",
      "\[IndentingNewLine]", 
     RowBox[{"GateChannel", "\[Rule]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"gateChannels", "[", 
         RowBox[{"[", 
          RowBox[{"Last", "[", 
           RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}], ",",
      "\[IndentingNewLine]", 
     RowBox[{"GateNoiseMemoryDepth", "\[Rule]", "1"}]}], 
    "\[IndentingNewLine]", "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"{", 
  RowBox[{
   RowBox[{"p", "[", "gs", "]"}], ",", 
   RowBox[{"t", "[", "gs", "]"}], ",", 
   RowBox[{"1", "-", 
    RowBox[{"F", "[", "gs", "]"}]}]}], "}"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"0.9997999995805747`", ",", "1.`", ",", "0.00010000020971268064`"}],
   "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"protocol", "=", 
   RowBox[{"RBProtocol", "[", 
    RowBox[{"mlist", ",", "5000", ",", "1", ",", 
     RowBox[{"TP", "[", "U", "]"}], ",", 
     RowBox[{"0.99", 
      RowBox[{"TP", "[", "U", "]"}]}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"survivals", "=", 
   RowBox[{"SimulateProtocol", "[", 
    RowBox[{"gs", ",", "protocol", ",", 
     RowBox[{
     "SimulationExportName", "\[Rule]", 
      "\"\<../data/overrotation_and_dephasing_model_survivals\>\""}]}], 
    "]"}]}], ";"}]}], "Input"],

Cell[BoxData["\<\"Saved ../data/overrotation_and_dephasing_model_survivals.h5\
\"\>"], "Print"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"Histogram", "[", 
    RowBox[{"#", ",", 
     RowBox[{"{", 
      RowBox[{"0", ",", "1", ",", "0.01"}], "}"}], ",", 
     RowBox[{"ImageSize", "\[Rule]", "400"}]}], "]"}], "&"}], "/@", 
  RowBox[{"survivals", "[", 
   RowBox[{"[", 
    RowBox[{"All", ",", "1", ",", "All", ",", "1"}], "]"}], "]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.98, 0}, {0.99, 4203},
        RoundingRadius->0], RectangleBox[{0.99, 0}, {1., 797},
        RoundingRadius->0]}, {}, {}}, {{}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.91, 0}, {0.92, 1},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 1},
        RoundingRadius->0], RectangleBox[{0.93, 0}, {0.9400000000000001, 6},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 34},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 106},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 370},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 1276},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 3206},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.86, 0}, {0.87, 3},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 3},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 1},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 4},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 28},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 24},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 64},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 138},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 238},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 424},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 806},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 1571},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 1696},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.71, 0}, {0.72, 3},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 1},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 1},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 2},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 3},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 4},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 4},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 7},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 4},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 9},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 15},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 24},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 20},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 29},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 52},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 64},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 87},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 100},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 146},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 190},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 243},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 329},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 359},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 512},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 656},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 890},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 1182},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 64},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.31, 0}, {0.32, 1},
        RoundingRadius->0], RectangleBox[{0.46, 0}, {0.47000000000000003, 1},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 1},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 1},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 1},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 1},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 2},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 4},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 1},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 4},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 2},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 7},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 4},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 6},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 6},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 6},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 11},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 11},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 18},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 13},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 21},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 24},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 32},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 35},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 48},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 47},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 62},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 65},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 61},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 83},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 84},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 85},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 127},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 133},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 168},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 187},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 236},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 227},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 316},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 331},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 364},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 420},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 492},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 546},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 627},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 78},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.2, 0}, {0.21, 1},
        RoundingRadius->0], RectangleBox[{0.24, 0}, {0.25, 1},
        RoundingRadius->0], RectangleBox[{0.27, 0}, {0.28, 4},
        RoundingRadius->0], RectangleBox[{0.3, 0}, {0.31, 1},
        RoundingRadius->0], RectangleBox[{0.31, 0}, {0.32, 1},
        RoundingRadius->0], RectangleBox[{0.32, 0}, {0.33, 5},
        RoundingRadius->0], RectangleBox[{0.34, 0}, {0.35000000000000003, 2},
        RoundingRadius->0], RectangleBox[{0.35000000000000003, 0}, {0.36, 4},
        RoundingRadius->0], RectangleBox[{0.36, 0}, {0.37, 2},
        RoundingRadius->0], RectangleBox[{0.37, 0}, {0.38, 4},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 4},
        RoundingRadius->0], RectangleBox[{0.39, 0}, {0.4, 2},
        RoundingRadius->0], RectangleBox[{0.4, 0}, {0.41000000000000003, 1},
        RoundingRadius->0], RectangleBox[{0.41000000000000003, 0}, {0.42, 6},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 3},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 5},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 10},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 8},
        RoundingRadius->0], RectangleBox[{0.46, 0}, {0.47000000000000003, 8},
        RoundingRadius->0], RectangleBox[{0.47000000000000003, 0}, {0.48, 8},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 8},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 8},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 14},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 12},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 19},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 11},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 19},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 15},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 15},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 21},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 25},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 23},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 22},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 31},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 35},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 33},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 33},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 31},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 42},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 55},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 49},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 59},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 64},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 71},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 69},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 75},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 104},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 91},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 99},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 99},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 100},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 120},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 147},
        RoundingRadius->0], 
       RectangleBox[{0.81, 0}, {0.8200000000000001, 139},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 145},
        RoundingRadius->0], 
       RectangleBox[{0.8300000000000001, 0}, {0.84, 166},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 171},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 181},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 220},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 183},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 232},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 245},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 247},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 269},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 319},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 309},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 365},
        RoundingRadius->0], 
       RectangleBox[{0.9500000000000001, 0}, {0.96, 110},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.08, 0}, {0.09, 1},
        RoundingRadius->0], RectangleBox[{0.09, 0}, {0.1, 4},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 2},
        RoundingRadius->0], RectangleBox[{0.11, 0}, {0.12, 5},
        RoundingRadius->0], RectangleBox[{0.12, 0}, {0.13, 8},
        RoundingRadius->0], RectangleBox[{0.13, 0}, {0.14, 5},
        RoundingRadius->0], RectangleBox[{0.14, 0}, {0.15, 9},
        RoundingRadius->0], RectangleBox[{0.15, 0}, {0.16, 4},
        RoundingRadius->0], RectangleBox[{0.16, 0}, {0.17, 9},
        RoundingRadius->0], RectangleBox[{0.17, 0}, {0.18, 10},
        RoundingRadius->0], RectangleBox[{0.18, 0}, {0.19, 7},
        RoundingRadius->0], RectangleBox[{0.19, 0}, {0.2, 18},
        RoundingRadius->0], RectangleBox[{0.2, 0}, {0.21, 12},
        RoundingRadius->0], RectangleBox[{0.21, 0}, {0.22, 13},
        RoundingRadius->0], RectangleBox[{0.22, 0}, {0.23, 14},
        RoundingRadius->0], RectangleBox[{0.23, 0}, {0.24, 12},
        RoundingRadius->0], RectangleBox[{0.24, 0}, {0.25, 10},
        RoundingRadius->0], RectangleBox[{0.25, 0}, {0.26, 17},
        RoundingRadius->0], RectangleBox[{0.26, 0}, {0.27, 16},
        RoundingRadius->0], RectangleBox[{0.27, 0}, {0.28, 19},
        RoundingRadius->0], RectangleBox[{0.28, 0}, {0.29, 17},
        RoundingRadius->0], RectangleBox[{0.29, 0}, {0.3, 21},
        RoundingRadius->0], RectangleBox[{0.3, 0}, {0.31, 13},
        RoundingRadius->0], RectangleBox[{0.31, 0}, {0.32, 15},
        RoundingRadius->0], RectangleBox[{0.32, 0}, {0.33, 17},
        RoundingRadius->0], RectangleBox[{0.33, 0}, {0.34, 23},
        RoundingRadius->0], 
       RectangleBox[{0.34, 0}, {0.35000000000000003, 26},
        RoundingRadius->0], 
       RectangleBox[{0.35000000000000003, 0}, {0.36, 23},
        RoundingRadius->0], RectangleBox[{0.36, 0}, {0.37, 26},
        RoundingRadius->0], RectangleBox[{0.37, 0}, {0.38, 36},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 19},
        RoundingRadius->0], RectangleBox[{0.39, 0}, {0.4, 38},
        RoundingRadius->0], RectangleBox[{0.4, 0}, {0.41000000000000003, 31},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 37},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 28},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 36},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 49},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 42},
        RoundingRadius->0], 
       RectangleBox[{0.46, 0}, {0.47000000000000003, 45},
        RoundingRadius->0], 
       RectangleBox[{0.47000000000000003, 0}, {0.48, 46},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 44},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 42},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 47},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 55},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 67},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 39},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 54},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 50},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 66},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 61},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 66},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 66},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 63},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 82},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 79},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 86},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 82},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 86},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 80},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 111},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 81},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 98},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 96},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 106},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 114},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 111},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 124},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 133},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 111},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 117},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 102},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 131},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 117},
        RoundingRadius->0], 
       RectangleBox[{0.81, 0}, {0.8200000000000001, 124},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 142},
        RoundingRadius->0], 
       RectangleBox[{0.8300000000000001, 0}, {0.84, 138},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 143},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 161},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 163},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 180},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 193},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 145},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 61},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.15, 0}, {0.16, 3},
        RoundingRadius->0], RectangleBox[{0.16, 0}, {0.17, 27},
        RoundingRadius->0], RectangleBox[{0.17, 0}, {0.18, 30},
        RoundingRadius->0], RectangleBox[{0.18, 0}, {0.19, 41},
        RoundingRadius->0], RectangleBox[{0.19, 0}, {0.2, 39},
        RoundingRadius->0], RectangleBox[{0.2, 0}, {0.21, 33},
        RoundingRadius->0], RectangleBox[{0.21, 0}, {0.22, 43},
        RoundingRadius->0], RectangleBox[{0.22, 0}, {0.23, 36},
        RoundingRadius->0], RectangleBox[{0.23, 0}, {0.24, 50},
        RoundingRadius->0], RectangleBox[{0.24, 0}, {0.25, 40},
        RoundingRadius->0], RectangleBox[{0.25, 0}, {0.26, 45},
        RoundingRadius->0], RectangleBox[{0.26, 0}, {0.27, 50},
        RoundingRadius->0], RectangleBox[{0.27, 0}, {0.28, 52},
        RoundingRadius->0], RectangleBox[{0.28, 0}, {0.29, 44},
        RoundingRadius->0], RectangleBox[{0.29, 0}, {0.3, 49},
        RoundingRadius->0], RectangleBox[{0.3, 0}, {0.31, 56},
        RoundingRadius->0], RectangleBox[{0.31, 0}, {0.32, 55},
        RoundingRadius->0], RectangleBox[{0.32, 0}, {0.33, 40},
        RoundingRadius->0], RectangleBox[{0.33, 0}, {0.34, 40},
        RoundingRadius->0], 
       RectangleBox[{0.34, 0}, {0.35000000000000003, 54},
        RoundingRadius->0], 
       RectangleBox[{0.35000000000000003, 0}, {0.36, 58},
        RoundingRadius->0], RectangleBox[{0.36, 0}, {0.37, 49},
        RoundingRadius->0], RectangleBox[{0.37, 0}, {0.38, 51},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 57},
        RoundingRadius->0], RectangleBox[{0.39, 0}, {0.4, 58},
        RoundingRadius->0], RectangleBox[{0.4, 0}, {0.41000000000000003, 59},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 66},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 66},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 63},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 72},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 65},
        RoundingRadius->0], 
       RectangleBox[{0.46, 0}, {0.47000000000000003, 75},
        RoundingRadius->0], 
       RectangleBox[{0.47000000000000003, 0}, {0.48, 76},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 81},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 78},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 72},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 77},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 59},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 73},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 88},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 69},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 98},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 97},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 80},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 84},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 98},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 82},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 91},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 79},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 85},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 95},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 97},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 89},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 96},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 113},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 91},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 125},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 97},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 123},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 99},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 99},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 113},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 131},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 121},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 98},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 118},
        RoundingRadius->0], 
       RectangleBox[{0.81, 0}, {0.8200000000000001, 120},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 126},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 16},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.26, 0}, {0.27, 26},
        RoundingRadius->0], RectangleBox[{0.27, 0}, {0.28, 90},
        RoundingRadius->0], RectangleBox[{0.28, 0}, {0.29, 92},
        RoundingRadius->0], RectangleBox[{0.29, 0}, {0.3, 109},
        RoundingRadius->0], RectangleBox[{0.3, 0}, {0.31, 104},
        RoundingRadius->0], RectangleBox[{0.31, 0}, {0.32, 102},
        RoundingRadius->0], RectangleBox[{0.32, 0}, {0.33, 100},
        RoundingRadius->0], RectangleBox[{0.33, 0}, {0.34, 104},
        RoundingRadius->0], 
       RectangleBox[{0.34, 0}, {0.35000000000000003, 97},
        RoundingRadius->0], 
       RectangleBox[{0.35000000000000003, 0}, {0.36, 93},
        RoundingRadius->0], RectangleBox[{0.36, 0}, {0.37, 102},
        RoundingRadius->0], RectangleBox[{0.37, 0}, {0.38, 98},
        RoundingRadius->0], RectangleBox[{0.38, 0}, {0.39, 96},
        RoundingRadius->0], RectangleBox[{0.39, 0}, {0.4, 94},
        RoundingRadius->0], 
       RectangleBox[{0.4, 0}, {0.41000000000000003, 100},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 90},
        RoundingRadius->0], RectangleBox[{0.42, 0}, {0.43, 111},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 112},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 104},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 83},
        RoundingRadius->0], 
       RectangleBox[{0.46, 0}, {0.47000000000000003, 102},
        RoundingRadius->0], 
       RectangleBox[{0.47000000000000003, 0}, {0.48, 108},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 109},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 123},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 108},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 114},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 122},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 124},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 114},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 99},
        RoundingRadius->0], 
       RectangleBox[{0.56, 0}, {0.5700000000000001, 118},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 99},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 126},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 114},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 119},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 124},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 136},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 110},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 125},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 124},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 119},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 126},
        RoundingRadius->0], 
       RectangleBox[{0.68, 0}, {0.6900000000000001, 122},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 117},
        RoundingRadius->0], 
       RectangleBox[{0.7000000000000001, 0}, {0.71, 120},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 118},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 53},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.42, 0}, {0.43, 253},
        RoundingRadius->0], RectangleBox[{0.43, 0}, {0.44, 328},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 343},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 350},
        RoundingRadius->0], 
       RectangleBox[{0.46, 0}, {0.47000000000000003, 332},
        RoundingRadius->0], 
       RectangleBox[{0.47000000000000003, 0}, {0.48, 388},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 356},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 350},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 353},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 335},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 337},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 348},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 342},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 352},
        RoundingRadius->0], 
       RectangleBox[{0.56, 0}, {0.5700000000000001, 233},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}]}], "}"}]], "Output"]
}, Open  ]]
}, Closed]],

Cell[CellGroupData[{

Cell["Pathological Error Model", "Section"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"\[Rho]r", "=", 
   RowBox[{
    RowBox[{"Unitary", "[", 
     RowBox[{"MatrixExp", "[", 
      RowBox[{
       RowBox[{"-", "I"}], " ", "0.1", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{"TP", "[", "X", "]"}], "+", 
          RowBox[{"TP", "[", "Y", "]"}]}], ")"}], "/", "2"}]}], "]"}], "]"}], 
    "[", 
    RowBox[{"TP", "[", "U", "]"}], "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"resetGate", "=", 
   RowBox[{"Super", "@", 
    RowBox[{"FunctionChannel", "[", 
     RowBox[{
      RowBox[{
       RowBox[{
        RowBox[{"Tr", "[", "#", "]"}], "\[Rho]r"}], "&"}], ",", 
      RowBox[{"InputDim", "\[Rule]", "2"}]}], "]"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"otherGate", "=", 
    RowBox[{"With", "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"p", "=", "0.01"}], "}"}], ",", 
      RowBox[{"Super", "@", 
       RowBox[{"FunctionChannel", "[", 
        RowBox[{
         RowBox[{
          RowBox[{
           RowBox[{"p", " ", 
            RowBox[{"Tr", "[", "#", "]"}], 
            RowBox[{
             RowBox[{"TP", "[", "I", "]"}], "/", "2"}]}], "+", 
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "-", "p"}], ")"}], " ", "#"}]}], "&"}], ",", 
         RowBox[{"InputDim", "\[Rule]", "2"}]}], "]"}]}]}], "]"}]}], ";"}], 
  "\[IndentingNewLine]", 
  RowBox[{"(*", 
   RowBox[{
    RowBox[{"otherGate", "=", 
     RowBox[{"With", "[", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"p", "=", "0.01"}], "}"}], ",", 
       RowBox[{"Super", "@", 
        RowBox[{"Kraus", "[", 
         RowBox[{"{", 
          RowBox[{
           RowBox[{
            RowBox[{"Sqrt", "[", 
             RowBox[{"1", "-", "p"}], "]"}], 
            RowBox[{"TP", "[", "I", "]"}]}], ",", 
           RowBox[{
            RowBox[{"Sqrt", "[", "p", "]"}], 
            RowBox[{"RandomUnitary", "[", "2", "]"}]}]}], "}"}], "]"}]}]}], 
      "]"}]}], ";"}], "*)"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"totalGate", "=", 
   RowBox[{"With", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"p", "=", "0.9"}], "}"}], ",", 
     RowBox[{
      RowBox[{"p", " ", "resetGate"}], "+", 
      RowBox[{
       RowBox[{"(", 
        RowBox[{"1", "-", "p"}], ")"}], " ", "otherGate"}]}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateNoise", "=", 
   RowBox[{"IndependentNoise", "[", "totalGate", "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"gateChannels", "=", 
    RowBox[{
     RowBox[{
      RowBox[{
       RowBox[{"Unitary", "[", 
        RowBox[{
         RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", "#",
          "]"}], "]"}], ".", 
       RowBox[{"First", "[", "gateNoise", "]"}]}], "&"}], "/@", 
     RowBox[{"Range", "[", "12", "]"}]}]}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"gs", "=", 
    RowBox[{"GateNoiseGateSet", "[", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{
      "GateSetName", "\[Rule]", "\"\<Multimodal survival distributions\>\""}],
       ",", "\[IndentingNewLine]", 
      RowBox[{"Size", "\[Rule]", "12"}], ",", "\[IndentingNewLine]", 
      RowBox[{"Dimension", "\[Rule]", "2"}], ",", "\[IndentingNewLine]", 
      RowBox[{"GateProduct", "\[Rule]", 
       RowBox[{"GateProduct", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
      "\[IndentingNewLine]", 
      RowBox[{"GateInverse", "\[Rule]", 
       RowBox[{"GateInverse", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
      "\[IndentingNewLine]", 
      RowBox[{"GateUnitary", "\[Rule]", 
       RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
      "\[IndentingNewLine]", 
      RowBox[{"GateNoise", "\[Rule]", "gateNoise"}], ",", 
      "\[IndentingNewLine]", 
      RowBox[{"GateChannel", "\[Rule]", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"gateChannels", "[", 
          RowBox[{"[", 
           RowBox[{"Last", "[", 
            RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}], 
      ",", "\[IndentingNewLine]", 
      RowBox[{"GateNoiseMemoryDepth", "\[Rule]", "1"}]}], 
     "\[IndentingNewLine]", "]"}]}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"protocol", "=", 
   RowBox[{"RBProtocol", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", "5", ",", "20", ",", "50", ",", "100"}], 
      "}"}], ",", "5000", ",", "1", ",", 
     RowBox[{"TP", "[", "U", "]"}], ",", 
     RowBox[{"0.99", 
      RowBox[{"TP", "[", "U", "]"}]}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"survivals", "=", 
   RowBox[{"SimulateProtocol", "[", 
    RowBox[{"gs", ",", "protocol", ",", 
     RowBox[{
     "SimulationExportName", "\[Rule]", 
      "\"\<../data/pathological_model_survivals\>\""}]}], "]"}]}], 
  ";"}]}], "Input"],

Cell[BoxData["\<\"Saved ../data/pathological_model_survivals.h5\"\>"], "Print"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"p", "[", "gs", "]"}], ",", 
   RowBox[{"t", "[", "gs", "]"}], ",", 
   RowBox[{"1", "-", 
    RowBox[{"F", "[", "gs", "]"}]}]}], "}"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0.09899999999999991`", ",", "0.9999999999999999`", ",", 
   "0.4505000000000001`"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{"Histogram", "[", 
    RowBox[{"#", ",", 
     RowBox[{"{", 
      RowBox[{"0", ",", "1", ",", "0.01"}], "}"}], ",", 
     RowBox[{"ImageSize", "\[Rule]", "400"}]}], "]"}], "&"}], "/@", 
  RowBox[{"survivals", "[", 
   RowBox[{"[", 
    RowBox[{"All", ",", "1", ",", "All", ",", "1"}], "]"}], 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"ListPlot", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"Flatten", "@", "#"}], ",", 
        RowBox[{"ConstantArray", "[", 
         RowBox[{"1", ",", 
          RowBox[{"Length", "@", 
           RowBox[{"Flatten", "@", "#"}]}]}], "]"}]}], "}"}], 
      "\[Transpose]"}], ",", 
     RowBox[{"Filling", "\[Rule]", "Axis"}], ",", 
     RowBox[{"ImageSize", "\[Rule]", "400"}], ",", 
     RowBox[{"PlotRange", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"0", ",", "1"}], "}"}], ",", 
        RowBox[{"{", 
         RowBox[{"0", ",", "1.5"}], "}"}]}], "}"}]}]}], "]"}], "&"}], "/@", 
  "survivals"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.1, 0}, {0.11, 824},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 1635},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 1697},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 844},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.01, 0}, {0.02, 139},
        RoundingRadius->0], RectangleBox[{0.05, 0}, {0.06, 275},
        RoundingRadius->0], RectangleBox[{0.06, 0}, {0.07, 308},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 134},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 275},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 1090},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 291},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 248},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 563},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 558},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 282},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 133},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 282},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 279},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 143},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0., 0}, {0.01, 67},
        RoundingRadius->0], RectangleBox[{0.01, 0}, {0.02, 74},
        RoundingRadius->0], RectangleBox[{0.04, 0}, {0.05, 182},
        RoundingRadius->0], RectangleBox[{0.05, 0}, {0.06, 328},
        RoundingRadius->0], RectangleBox[{0.06, 0}, {0.07, 46},
        RoundingRadius->0], RectangleBox[{0.09, 0}, {0.1, 110},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 20},
        RoundingRadius->0], 
       RectangleBox[{0.4, 0}, {0.41000000000000003, 236},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 49},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 455},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 667},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 35},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 469},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 36},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 676},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 440},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 51},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 244},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 26},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 114},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 44},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 311},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 175},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 80},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 65},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0., 0}, {0.01, 62},
        RoundingRadius->0], RectangleBox[{0.01, 0}, {0.02, 73},
        RoundingRadius->0], RectangleBox[{0.04, 0}, {0.05, 208},
        RoundingRadius->0], RectangleBox[{0.05, 0}, {0.06, 332},
        RoundingRadius->0], RectangleBox[{0.06, 0}, {0.07, 39},
        RoundingRadius->0], RectangleBox[{0.09, 0}, {0.1, 122},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 15},
        RoundingRadius->0], 
       RectangleBox[{0.4, 0}, {0.41000000000000003, 231},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 58},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 452},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 650},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 21},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 504},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 46},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 647},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 453},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 33},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 213},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 15},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 97},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 48},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 332},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 194},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 87},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 68},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0., 0}, {0.01, 58},
        RoundingRadius->0], RectangleBox[{0.01, 0}, {0.02, 75},
        RoundingRadius->0], RectangleBox[{0.04, 0}, {0.05, 171},
        RoundingRadius->0], RectangleBox[{0.05, 0}, {0.06, 322},
        RoundingRadius->0], RectangleBox[{0.06, 0}, {0.07, 46},
        RoundingRadius->0], RectangleBox[{0.09, 0}, {0.1, 117},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 26},
        RoundingRadius->0], 
       RectangleBox[{0.4, 0}, {0.41000000000000003, 237},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 44},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 436},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 672},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 31},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 506},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 39},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 629},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 489},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 35},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 234},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 27},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 116},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 41},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 304},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 194},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 89},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 62},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0., 0}, {0.01, 52},
        RoundingRadius->0], RectangleBox[{0.01, 0}, {0.02, 86},
        RoundingRadius->0], RectangleBox[{0.04, 0}, {0.05, 214},
        RoundingRadius->0], RectangleBox[{0.05, 0}, {0.06, 277},
        RoundingRadius->0], RectangleBox[{0.06, 0}, {0.07, 42},
        RoundingRadius->0], RectangleBox[{0.09, 0}, {0.1, 123},
        RoundingRadius->0], RectangleBox[{0.1, 0}, {0.11, 26},
        RoundingRadius->0], 
       RectangleBox[{0.4, 0}, {0.41000000000000003, 220},
        RoundingRadius->0], 
       RectangleBox[{0.41000000000000003, 0}, {0.42, 43},
        RoundingRadius->0], RectangleBox[{0.44, 0}, {0.45, 446},
        RoundingRadius->0], RectangleBox[{0.45, 0}, {0.46, 659},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 45},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 493},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 37},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 634},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 488},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 43},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 207},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 23},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 111},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 48},
        RoundingRadius->0], 
       RectangleBox[{0.93, 0}, {0.9400000000000001, 344},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 182},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 80},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 77},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}]}], "}"}]], "Output"]
}, Open  ]]
}, Closed]],

Cell[CellGroupData[{

Cell["X/2 Only", "Section"],

Cell["\<\
These are the necessary amplitudes (in MHz) for each of 10 1ns pulses steps \
to get a perfect \[Pi]/2 rotation about x:\
\>", "Text"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"rawPulse", "=", 
  RowBox[{"GaussianTailsPulse", "[", 
   RowBox[{"0.001", ",", "0.01", ",", "0.005", ",", 
    RowBox[{"Area", "\[Rule]", 
     RowBox[{"1", "/", "4"}]}]}], "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"X2unitary", "=", 
   RowBox[{
    RowBox[{"Last", "@", 
     RowBox[{"Unitaries", "@", 
      RowBox[{"PulseSim", "[", 
       RowBox[{
        RowBox[{"2", "\[Pi]", " ", "0", "*", 
         RowBox[{
          RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], ",", 
        RowBox[{"{", 
         RowBox[{"rawPulse", ",", 
          RowBox[{"{", 
           RowBox[{
            RowBox[{"2", "\[Pi]", " ", 
             RowBox[{
              RowBox[{"TP", "[", "X", "]"}], "/", "2"}]}], ",", 
            RowBox[{"2", "\[Pi]", " ", 
             RowBox[{
              RowBox[{"TP", "[", "Y", "]"}], "/", "2"}]}]}], "}"}]}], "}"}]}],
        "]"}]}]}], "//", "Chop"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"X2unitary", "//", "MatrixForm"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "0.25871966577139427`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "1.5651627700276785`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "6.606075956374566`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "19.452776781521813`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "39.964430659200744`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "57.282193999989516`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "57.282193999989516`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "39.96443065920073`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "19.45277678152182`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "6.606075956374562`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.001`", ",", "1.5651627700276807`", ",", "0.`"}], "}"}]}], 
  "}"}]], "Output"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"0.7071067811865474`", 
      RowBox[{"0.`", "\[VeryThinSpace]", "-", 
       RowBox[{"0.7071067811865479`", " ", "\[ImaginaryI]"}]}]},
     {
      RowBox[{"0.`", "\[VeryThinSpace]", "-", 
       RowBox[{"0.707106781186548`", " ", "\[ImaginaryI]"}]}], 
      "0.7071067811865477`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[BoxData[{
 RowBox[{
  RowBox[{"ClearAll", "@", "X2channel"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"X2channel", "[", 
   RowBox[{
   "\[Epsilon]mean_", ",", "\[Epsilon]std_", ",", "\[Delta]\[Omega]mean_", 
    ",", "\[Delta]\[Omega]std_", ",", "T2_"}], "]"}], ":=", 
  RowBox[{
   RowBox[{"X2channel", "[", 
    RowBox[{
    "\[Epsilon]mean", ",", "\[Epsilon]std", ",", "\[Delta]\[Omega]mean", ",", 
     "\[Delta]\[Omega]std"}], "]"}], "=", 
   RowBox[{"Module", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Delta]\[Omega]", ",", "\[Epsilon]"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"Last", "@", 
      RowBox[{"Superoperators", "@", 
       RowBox[{"PulseSim", "[", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"LindbladForm", "[", 
          RowBox[{
           RowBox[{"2", "\[Pi]", " ", "\[Delta]\[Omega]", "*", 
            RowBox[{
             RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], ",", 
           RowBox[{"{", 
            RowBox[{
             RowBox[{"1", "/", 
              RowBox[{"Sqrt", "[", "T2", "]"}]}], 
             RowBox[{
              RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], "}"}]}], "]"}], 
         ",", "\[IndentingNewLine]", 
         RowBox[{"{", 
          RowBox[{"rawPulse", ",", 
           RowBox[{"\[Epsilon]", 
            RowBox[{"{", 
             RowBox[{
              RowBox[{"2", "\[Pi]", " ", 
               RowBox[{
                RowBox[{"TP", "[", "X", "]"}], "/", "2"}]}], ",", 
              RowBox[{"2", "\[Pi]", " ", 
               RowBox[{
                RowBox[{"TP", "[", "Y", "]"}], "/", "2"}]}]}], "}"}]}]}], 
          "}"}], ",", "\[IndentingNewLine]", 
         StyleBox[
          RowBox[{"RandomMultinormalParameterDistribution", "[", 
           RowBox[{
            RowBox[{"{", 
             RowBox[{"\[Epsilon]mean", ",", "\[Delta]\[Omega]mean"}], "}"}], 
            ",", 
            RowBox[{"DiagonalMatrix", "[", 
             RowBox[{
              RowBox[{"{", 
               RowBox[{"\[Epsilon]std", ",", "\[Delta]\[Omega]std"}], "}"}], 
              "^", "2"}], "]"}], ",", 
            RowBox[{"{", 
             RowBox[{"\[Epsilon]", ",", "\[Delta]\[Omega]"}], "}"}], ",", 
            "100"}], "]"}], "Input"]}], "\[IndentingNewLine]", "]"}]}]}]}], 
    "\[IndentingNewLine]", "]"}]}]}]}], "Input"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"X2choice", "=", 
   RowBox[{"X2channel", "[", 
    RowBox[{"1.001", ",", "0.01", ",", "0.001", ",", "0.01", ",", "100"}], 
    "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"Re", "@", 
  RowBox[{"AverageGateFidelity", "[", 
   RowBox[{
    RowBox[{"Unitary", "[", "X2choice", "]"}], ",", "example"}], 
   "]"}]}]}], "Input"],

Cell[BoxData[
 RowBox[{"Re", "[", 
  RowBox[{"AverageGateFidelity", "[", 
   RowBox[{
    RowBox[{"\<\"Super\"\>", "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{
          RowBox[{"0.4989376913635425`", "\[VeryThinSpace]", "+", 
           RowBox[{"9.11815793152317`*^-23", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{"0.00003519691841870268`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.49991309758398`", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{"0.000035196918418702675`", "\[VeryThinSpace]", "+", 
           RowBox[{"0.49991309758397995`", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{"0.5010623086364577`", "\[VeryThinSpace]", "+", 
           RowBox[{"4.270429835596486`*^-22", " ", "\[ImaginaryI]"}]}]}], 
         "}"}], ",", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{"0.00003004259293401474`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.4999155946360185`", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{"0.49890121679142074`", "\[VeryThinSpace]", "+", 
           RowBox[{"0.00006576322424396107`", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{"0.5010436338685813`", "\[VeryThinSpace]", "-", 
           RowBox[{"5.113054689865564`*^-6", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{
           RowBox[{"-", "0.000030042592934014755`"}], "+", 
           RowBox[{"0.4999155946360185`", " ", "\[ImaginaryI]"}]}]}], "}"}], 
        ",", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{"0.00003004259293401474`", "\[VeryThinSpace]", "+", 
           RowBox[{"0.4999155946360186`", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{"0.5010436338685813`", "\[VeryThinSpace]", "+", 
           RowBox[{"5.1130546898655665`*^-6", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{"0.49890121679142085`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.00006576322424396107`", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{
           RowBox[{"-", "0.000030042592934014755`"}], "-", 
           RowBox[{"0.49991559463601865`", " ", "\[ImaginaryI]"}]}]}], "}"}], 
        ",", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{"0.5010623086364577`", "\[VeryThinSpace]", "-", 
           RowBox[{"3.3294853540816546`*^-22", " ", "\[ImaginaryI]"}]}], ",", 
          
          RowBox[{
           RowBox[{"-", "0.000035196918418702695`"}], "+", 
           RowBox[{"0.49991309758398`", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{
           RowBox[{"-", "0.00003519691841870269`"}], "-", 
           RowBox[{"0.49991309758397984`", " ", "\[ImaginaryI]"}]}], ",", 
          RowBox[{"0.49893769136354266`", "\[VeryThinSpace]", "+", 
           RowBox[{"1.6124557385860905`*^-22", " ", "\[ImaginaryI]"}]}]}], 
         "}"}]}], "}"}], ",", "\<\"<params>\"\>"}], "]"}], ",", "example"}], 
   "]"}], "]"}]], "Output"]
}, Open  ]],

Cell[BoxData[{
 RowBox[{
  RowBox[{"group", "=", 
   RowBox[{"Dot", "@@@", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"Reverse", "/@", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"W", "[", "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"X2", ",", "X2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"Z2", ",", "X2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"Z2", ",", "X2", ",", "Z"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"ZM2", ",", "X2", ",", "Z"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"Z", ",", "X2", ",", "Z2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"ZM2", ",", "X2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"X2", ",", "ZM2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"Z", ",", "X2", ",", "ZM2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"X2", ",", "Z2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"X2", ",", "X2", ",", "Z"}], "]"}], ",", 
         RowBox[{"W", "[", "Z", "]"}]}], "}"}]}], "/.", 
      RowBox[{
       RowBox[{"W", "[", "]"}], "\[Rule]", "Id"}]}], ")"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateChannels", "=", 
   RowBox[{"Super", "/@", 
    RowBox[{"(", 
     RowBox[{"group", "/.", 
      RowBox[{"{", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"X2", "->", 
         RowBox[{"X2channel", "[", 
          RowBox[{
          "1.001", ",", "0.01", ",", "0.001", ",", "0.01", ",", "100"}], 
          "]"}]}], ",", "\[IndentingNewLine]", 
        RowBox[{"Z2", "\[Rule]", 
         RowBox[{"Unitary", "[", 
          RowBox[{"MatrixExp", "[", 
           RowBox[{
            RowBox[{"-", "I"}], " ", 
            RowBox[{"(", 
             RowBox[{"\[Pi]", "/", "2"}], ")"}], 
            RowBox[{
             RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], "]"}], "]"}]}], ",",
         "\[IndentingNewLine]", 
        RowBox[{"Z", "\[Rule]", 
         RowBox[{"Unitary", "[", 
          RowBox[{"MatrixExp", "[", 
           RowBox[{
            RowBox[{"-", "I"}], " ", 
            RowBox[{"(", "\[Pi]", ")"}], 
            RowBox[{
             RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], "]"}], "]"}]}], ",",
         "\[IndentingNewLine]", 
        RowBox[{"ZM2", "\[Rule]", 
         RowBox[{"Unitary", "[", 
          RowBox[{"MatrixExp", "[", 
           RowBox[{"I", " ", 
            RowBox[{"(", 
             RowBox[{"\[Pi]", "/", "2"}], ")"}], 
            RowBox[{
             RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], "]"}], "]"}]}], ",",
         "\[IndentingNewLine]", 
        RowBox[{"Id", "\[Rule]", 
         RowBox[{"Unitary", "[", 
          RowBox[{"IdentityMatrix", "[", "2", "]"}], "]"}]}]}], 
       "\[IndentingNewLine]", "}"}]}], ")"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateUnitaries", "=", 
   RowBox[{"First", "/@", 
    RowBox[{"(", 
     RowBox[{"group", "/.", 
      RowBox[{"{", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"X2", "->", 
         RowBox[{"Unitary", "[", 
          RowBox[{"MatrixExp", "[", 
           RowBox[{
            RowBox[{"-", "I"}], " ", 
            RowBox[{"(", 
             RowBox[{"\[Pi]", "/", "2"}], ")"}], 
            RowBox[{
             RowBox[{"TP", "[", "X", "]"}], "/", "2"}]}], "]"}], "]"}]}], ",",
         "\[IndentingNewLine]", 
        RowBox[{"Z2", "\[Rule]", 
         RowBox[{"Unitary", "[", 
          RowBox[{"MatrixExp", "[", 
           RowBox[{
            RowBox[{"-", "I"}], " ", 
            RowBox[{"(", 
             RowBox[{"\[Pi]", "/", "2"}], ")"}], 
            RowBox[{
             RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], "]"}], "]"}]}], ",",
         "\[IndentingNewLine]", 
        RowBox[{"Z", "\[Rule]", 
         RowBox[{"Unitary", "[", 
          RowBox[{"MatrixExp", "[", 
           RowBox[{
            RowBox[{"-", "I"}], " ", 
            RowBox[{"(", "\[Pi]", ")"}], 
            RowBox[{
             RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], "]"}], "]"}]}], ",",
         "\[IndentingNewLine]", 
        RowBox[{"ZM2", "\[Rule]", 
         RowBox[{"Unitary", "[", 
          RowBox[{"MatrixExp", "[", 
           RowBox[{"I", " ", 
            RowBox[{"(", 
             RowBox[{"\[Pi]", "/", "2"}], ")"}], 
            RowBox[{
             RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], "]"}], "]"}]}], ",",
         "\[IndentingNewLine]", 
        RowBox[{"Id", "\[Rule]", 
         RowBox[{"Unitary", "[", 
          RowBox[{"IdentityMatrix", "[", "2", "]"}], "]"}]}]}], 
       "\[IndentingNewLine]", "}"}]}], ")"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateNoises", "=", 
   RowBox[{
    RowBox[{
     RowBox[{"Super", "[", 
      RowBox[{
       RowBox[{"Unitary", "[", 
        RowBox[{"#2", "\[ConjugateTranspose]"}], "]"}], ".", "#1"}], "]"}], 
     "&"}], "@@@", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"gateChannels", ",", "gateUnitaries"}], "}"}], 
      "\[Transpose]"}], ")"}]}]}], ";"}]}], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"GateProduct", "[", "$gaussianQubitGateSet", "]"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
      "1", ",", "2", ",", "3", ",", "4", ",", "5", ",", "6", ",", "7", ",", 
       "8", ",", "9", ",", "10", ",", "11", ",", "12"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "2", ",", "4", ",", "6", ",", "1", ",", "8", ",", "9", ",", "11", ",", 
       "12", ",", "3", ",", "7", ",", "10", ",", "5"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "3", ",", "5", ",", "1", ",", "7", ",", "2", ",", "10", ",", "4", ",", 
       "9", ",", "8", ",", "6", ",", "12", ",", "11"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "4", ",", "1", ",", "9", ",", "2", ",", "12", ",", "3", ",", "10", ",", 
       "5", ",", "6", ",", "11", ",", "7", ",", "8"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "5", ",", "7", ",", "10", ",", "3", ",", "9", ",", "8", ",", "12", ",", 
       "11", ",", "1", ",", "4", ",", "6", ",", "2"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "6", ",", "8", ",", "2", ",", "11", ",", "4", ",", "7", ",", "1", ",", 
       "3", ",", "12", ",", "9", ",", "5", ",", "10"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "7", ",", "3", ",", "8", ",", "5", ",", "11", ",", "1", ",", "6", ",", 
       "2", ",", "10", ",", "12", ",", "4", ",", "9"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "8", ",", "11", ",", "7", ",", "6", ",", "3", ",", "12", ",", "5", ",", 
       "10", ",", "2", ",", "1", ",", "9", ",", "4"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "9", ",", "12", ",", "4", ",", "10", ",", "1", ",", "11", ",", "2", ",",
        "6", ",", "5", ",", "3", ",", "8", ",", "7"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "10", ",", "9", ",", "5", ",", "12", ",", "7", ",", "4", ",", "3", ",", 
       "1", ",", "11", ",", "8", ",", "2", ",", "6"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "11", ",", "6", ",", "12", ",", "8", ",", "10", ",", "2", ",", "9", ",",
        "4", ",", "7", ",", "5", ",", "1", ",", "3"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "12", ",", "10", ",", "11", ",", "9", ",", "6", ",", "5", ",", "8", ",",
        "7", ",", "4", ",", "2", ",", "3", ",", "1"}], "}"}]}], "}"}], 
   "\[LeftDoubleBracket]", 
   RowBox[{"#1", ",", "#2"}], "\[RightDoubleBracket]"}], "&"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"?", "GateChannel"}]], "Input"],

Cell[BoxData[
 StyleBox["\<\"\!\(\*TagBox[StyleBox[\\\"GateChannel\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"]], DisplayForm]\) is a \
\!\(\*TagBox[StyleBox[\\\"GateNoiseGateSet\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"]], DisplayForm]\) key storing a function  \!\
\(\*TagBox[StyleBox[\\\"f[\\\", \\\"Input\\\", Rule[FontFamily, \\\"Courier\\\
\"]], DisplayForm]\)\!\(\*TagBox[StyleBox[\\\"i1\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"], Rule[FontWeight, \\\"Plain\\\"], \
Rule[FontSlant, \\\"Italic\\\"]], DisplayForm]\)\!\(\*TagBox[StyleBox[\\\",\\\
\", \\\"Input\\\", Rule[FontFamily, \\\"Courier\\\"]], \
DisplayForm]\)\!\(\*TagBox[StyleBox[\\\"i2\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"], Rule[FontWeight, \\\"Plain\\\"], \
Rule[FontSlant, \\\"Italic\\\"]], DisplayForm]\)\!\(\*TagBox[StyleBox[\\\",\\\
\", \\\"Input\\\", Rule[FontFamily, \\\"Courier\\\"]], \
DisplayForm]\)\!\(\*TagBox[StyleBox[\\\"...\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"], Rule[FontWeight, \\\"Plain\\\"], \
Rule[FontSlant, \\\"Italic\\\"]], DisplayForm]\)\!\(\*TagBox[StyleBox[\\\",\\\
\", \\\"Input\\\", Rule[FontFamily, \\\"Courier\\\"]], \
DisplayForm]\)\!\(\*TagBox[StyleBox[\\\"in\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"], Rule[FontWeight, \\\"Plain\\\"], \
Rule[FontSlant, \\\"Italic\\\"]], DisplayForm]\)\!\(\*TagBox[StyleBox[\\\"]\\\
\", \\\"Input\\\", Rule[FontFamily, \\\"Courier\\\"]], DisplayForm]\) that \
returns an object satisfying \!\(\*TagBox[StyleBox[\\\"ChannelPulseQ\\\", \
\\\"Input\\\", Rule[FontFamily, \\\"Courier\\\"]], DisplayForm]\) given the \
gate index history \!\(\*TagBox[StyleBox[\\\"i1,i2,...,in\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"], Rule[FontWeight, \\\"Plain\\\"], \
Rule[FontSlant, \\\"Italic\\\"]], DisplayForm]\) where \!\(\*TagBox[StyleBox[\
\\\"in\\\", \\\"Input\\\", Rule[FontFamily, \\\"Courier\\\"], \
Rule[FontWeight, \\\"Plain\\\"], Rule[FontSlant, \\\"Italic\\\"]], \
DisplayForm]\) is the most recent. This channel is the composition of the \
ideal gate \!\(\*TagBox[StyleBox[\\\"GateUnitary\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"]], DisplayForm]\) with index \
\!\(\*TagBox[StyleBox[\\\"in\\\", \\\"Input\\\", Rule[FontFamily, \\\"Courier\
\\\"], Rule[FontWeight, \\\"Plain\\\"], Rule[FontSlant, \\\"Italic\\\"]], \
DisplayForm]\) followed by the \!\(\*TagBox[StyleBox[\\\"GateNoise\\\", \
\\\"Input\\\", Rule[FontFamily, \\\"Courier\\\"]], DisplayForm]\) error, \
where the \!\(\*TagBox[StyleBox[\\\"GateNoise\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"]], DisplayForm]\) is allowed to depend on \
the indices\!\(\*TagBox[StyleBox[\\\" i1,i2,...,i(n-1)\\\", \\\"Input\\\", \
Rule[FontFamily, \\\"Courier\\\"], Rule[FontWeight, \\\"Plain\\\"], \
Rule[FontSlant, \\\"Italic\\\"]], DisplayForm]\) that preceded the current \
gate.\"\>", "MSG"]], "Print", "PrintUsage",
 CellTags->"Info993710147481-3174582"]
}, Open  ]],

Cell[BoxData["TestGateSetProducts"], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"gs", "=", 
   RowBox[{"GateNoiseGateSet", "[", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"GateSetName", "\[Rule]", "\"\<Depolarizing Noise\>\""}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"Size", "\[Rule]", "12"}], ",", "\[IndentingNewLine]", 
     RowBox[{"Dimension", "\[Rule]", "2"}], ",", "\[IndentingNewLine]", 
     RowBox[{"GateProduct", "\[Rule]", 
      RowBox[{"GateProduct", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateInverse", "\[Rule]", 
      RowBox[{"GateInverse", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateUnitary", "\[Rule]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"gateUnitaries", "[", 
         RowBox[{"[", 
          RowBox[{"Last", "[", 
           RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}], ",",
      "\[IndentingNewLine]", 
     RowBox[{"GateNoise", "\[Rule]", "gateNoise"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateChannel", "\[Rule]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"gateChannels", "[", 
         RowBox[{"[", 
          RowBox[{"Last", "[", 
           RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}]}], 
    "\[IndentingNewLine]", "]"}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"TestGateSetProducts", "[", "gs", "]"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{
    "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1",
      ",", "1", ",", "1", ",", "1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "0", ",", "0", ",", 
     FractionBox["1", "4"]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "0", ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "0", ",", "0", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "0", ",", "0", ",", 
     FractionBox["1", "4"], ",", 
     RowBox[{
      FractionBox["1", "4"], " ", 
      SuperscriptBox[
       RowBox[{"Abs", "[", 
        RowBox[{
         RowBox[{"-", 
          FractionBox[
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "-", "\[ImaginaryI]"}], ")"}], " ", 
            SuperscriptBox["\[ExponentialE]", 
             RowBox[{"-", 
              FractionBox[
               RowBox[{"\[ImaginaryI]", " ", "\[Pi]"}], "4"]}]]}], 
           SqrtBox["2"]]}], "-", 
         FractionBox[
          RowBox[{
           RowBox[{"(", 
            RowBox[{"1", "+", "\[ImaginaryI]"}], ")"}], " ", 
           SuperscriptBox["\[ExponentialE]", 
            FractionBox[
             RowBox[{"\[ImaginaryI]", " ", "\[Pi]"}], "4"]]}], 
          SqrtBox["2"]]}], "]"}], "2"]}], ",", "0", ",", 
     RowBox[{
      FractionBox["1", "4"], " ", 
      SuperscriptBox[
       RowBox[{"Abs", "[", 
        RowBox[{
         RowBox[{"-", 
          FractionBox[
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "-", "\[ImaginaryI]"}], ")"}], " ", 
            SuperscriptBox["\[ExponentialE]", 
             RowBox[{"-", 
              FractionBox[
               RowBox[{"\[ImaginaryI]", " ", "\[Pi]"}], "4"]}]]}], 
           SqrtBox["2"]]}], "-", 
         FractionBox[
          RowBox[{
           RowBox[{"(", 
            RowBox[{"1", "+", "\[ImaginaryI]"}], ")"}], " ", 
           SuperscriptBox["\[ExponentialE]", 
            FractionBox[
             RowBox[{"\[ImaginaryI]", " ", "\[Pi]"}], "4"]]}], 
          SqrtBox["2"]]}], "]"}], "2"]}], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", 
     RowBox[{
      FractionBox["1", "4"], " ", 
      SuperscriptBox[
       RowBox[{"Abs", "[", 
        RowBox[{
         FractionBox[
          RowBox[{
           RowBox[{"(", 
            RowBox[{"1", "-", "\[ImaginaryI]"}], ")"}], " ", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"-", 
             FractionBox[
              RowBox[{"\[ImaginaryI]", " ", "\[Pi]"}], "4"]}]]}], 
          SqrtBox["2"]], "+", 
         FractionBox[
          RowBox[{
           RowBox[{"(", 
            RowBox[{"1", "+", "\[ImaginaryI]"}], ")"}], " ", 
           SuperscriptBox["\[ExponentialE]", 
            FractionBox[
             RowBox[{"\[ImaginaryI]", " ", "\[Pi]"}], "4"]]}], 
          SqrtBox["2"]]}], "]"}], "2"]}], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     RowBox[{
      FractionBox["1", "4"], " ", 
      SuperscriptBox[
       RowBox[{"Abs", "[", 
        RowBox[{
         RowBox[{"-", 
          FractionBox[
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "-", "\[ImaginaryI]"}], ")"}], " ", 
            SuperscriptBox["\[ExponentialE]", 
             RowBox[{"-", 
              FractionBox[
               RowBox[{"\[ImaginaryI]", " ", "\[Pi]"}], "4"]}]]}], 
           SqrtBox["2"]]}], "-", 
         FractionBox[
          RowBox[{
           RowBox[{"(", 
            RowBox[{"1", "+", "\[ImaginaryI]"}], ")"}], " ", 
           SuperscriptBox["\[ExponentialE]", 
            FractionBox[
             RowBox[{"\[ImaginaryI]", " ", "\[Pi]"}], "4"]]}], 
          SqrtBox["2"]]}], "]"}], "2"]}], ",", "1", ",", 
     FractionBox["1", "4"], ",", "0", ",", "0", ",", 
     FractionBox["1", "4"]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     FractionBox["1", "4"]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     FractionBox["1", "4"], ",", "0", ",", "1", ",", "0", ",", 
     FractionBox["1", "4"]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     FractionBox["1", "4"], ",", "0", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "1", ",", 
     FractionBox["1", "4"]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", 
     FractionBox["1", "4"], ",", "1"}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"nm", "=", 
   RowBox[{"NoiseModel", "[", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"DistortionOperator", "\[Rule]", 
      RowBox[{"FastExponentialDistortion", "[", 
       RowBox[{"0.003", ",", "10", ",", "0"}], "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"DistortionMultiplier", "\[Rule]", "10"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateNoise", "\[Rule]", 
      RowBox[{"IndependentNoise", "[", 
       RowBox[{"Kraus", "[", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           RowBox[{"Sqrt", "[", "0.9999", "]"}], " ", 
           RowBox[{"TP", "[", "I", "]"}]}], ",", 
          RowBox[{
           RowBox[{"Sqrt", "[", "0.0001", "]"}], " ", 
           RowBox[{"TP", "[", "Z", "]"}]}]}], "}"}], "]"}], "]"}]}]}], 
    "\[IndentingNewLine]", "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"gs", "=", 
  RowBox[{"CompileGateNoise", "[", 
   RowBox[{"$gaussianQubitGateSet", ",", "nm", ",", "2"}], "]"}]}]}], "Input"],

Cell[BoxData[
 RowBox[{"\<\"GateNoiseGateSet\"\>", "[", 
  RowBox[{
   RowBox[{
   "GateSetName", "\[Rule]", "\<\"\\\"Gaussian 2-design on Qubits\\\"\"\>"}], 
   ",", 
   RowBox[{"Dimension", "\[Rule]", "2"}], ",", 
   RowBox[{"Size", "\[Rule]", "12"}], ",", "\<\"...\"\>"}], "]"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"protocol", "=", 
   RowBox[{"RBProtocol", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
      "1", ",", "10", ",", "50", ",", "100", ",", "500", ",", "1000", ",", 
       "1500", ",", "2000", ",", "3000", ",", "5000", ",", "10000", ",", 
       "20000"}], "}"}], ",", "100", ",", "1", ",", 
     RowBox[{"TP", "[", "U", "]"}], ",", 
     RowBox[{"TP", "[", "U", "]"}]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"survivals", "=", 
   RowBox[{"SimulateProtocol", "[", 
    RowBox[{"gs", ",", "protocol", ",", 
     RowBox[{
     "SimulationExportName", "\[Rule]", "\"\<../data/high_fidelity\>\""}]}], 
    "]"}]}], ";"}]}], "Input"],

Cell[BoxData["\<\"Saved ../data/high_fidelity.h5\"\>"], "Print"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"model", "=", 
   RowBox[{
    RowBox[{
     RowBox[{"(", 
      RowBox[{"A", "-", "B"}], ")"}], 
     SuperscriptBox["p", "m"]}], "+", "B"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"normalizedSurvivals", "=", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"SequenceLengths", "[", "protocol", "]"}], ",", 
      RowBox[{
       RowBox[{"Total", "[", 
        RowBox[{
         RowBox[{"survivals", "[", 
          RowBox[{"[", 
           RowBox[{"All", ",", "1", ",", "All", ",", "1"}], "]"}], "]"}], ",", 
         RowBox[{"{", "2", "}"}]}], "]"}], "/", "100"}]}], "}"}], 
    "\[Transpose]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"fit", "=", 
   RowBox[{"FindFit", "[", 
    RowBox[{"normalizedSurvivals", ",", 
     RowBox[{"{", 
      RowBox[{"model", ",", 
       RowBox[{"0", "<", "p", "<", "1"}], ",", 
       RowBox[{"0", "<", "A", "<", "1"}], ",", 
       RowBox[{"0", "<", "B", "<", "1"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"A", ",", "B", ",", "p"}], "}"}], ",", "m"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"modelf", "=", 
   RowBox[{"Function", "[", 
    RowBox[{
     RowBox[{"{", "m", "}"}], ",", 
     RowBox[{"Evaluate", "[", 
      RowBox[{"model", "/.", "fit"}], "]"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Print", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"A", ",", "B", ",", "p"}], "}"}], "/.", "fit"}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Plot", "[", 
  RowBox[{
   RowBox[{"modelf", "[", "m", "]"}], ",", 
   RowBox[{"{", 
    RowBox[{"m", ",", "0", ",", "20000"}], "}"}], ",", 
   RowBox[{"Epilog", "\[Rule]", 
    RowBox[{"Map", "[", 
     RowBox[{"Point", ",", "normalizedSurvivals"}], "]"}]}], ",", 
   RowBox[{"PlotRange", "\[Rule]", 
    RowBox[{"{", 
     RowBox[{"All", ",", 
      RowBox[{"{", 
       RowBox[{"0", ",", "1"}], "}"}]}], "}"}]}]}], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0.982911523186985`", ",", "0.4999738862680617`", ",", 
   "0.9998428930716474`"}], "}"}]], "Print"],

Cell[BoxData[
 GraphicsBox[{{{}, {}, 
    TagBox[
     {RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6], Opacity[
      1.], LineBox[CompressedData["
1:eJwVz2s0lAkcBvBh10Gii5RyS1ebpdKiUP9/aUqRlRW6UC5NOU3SoVaYHCV6
zVGRRG6TLTuTNqRUWiSDhGjellyqd5gZzZhLjqNcYvfdD895zu/D8+GxDj3l
e1SbwWDsofN/l6x1o7rrXCC4aXO7dowGOowC8QZhii/0rO03RGkgwJv/4jBh
g37qvMRwtgbMIjsZTGIj3gs+MdoSroFpjoTnTHjgxOpvi3j7NKC+pP6yjghE
XtsK63POGpARozl2xHHMCpDk9k+pYcDgq9KWiEVlr2fyoVQ1ZJwJ9/bPIPCY
aZpPY7YKlAFORq4zmThSaMlZL1SCkaLrbv6CXEx6L/xmpa0EYq59j3FZEVbf
j44rcBuGn71fhwUpi1E0SkzWXVUAfHzPd3UvwR2jXs2WlBxEfMemC2kCbMvK
K27dIoc5q5bUh+TeR8tI4dy8is8Qq+ohPBVleDZB80hr2WdorJybw0qqQPue
0kLfO0NgdOVgd59OJW4vMGbHrRuCk4ZZHP/SR5jXRxlerpaBNeva/IJfqpBv
1pez5zcZDIfFGck/PsGpofrnIQopjJ3nno1IeoZutS6sGUIKftFH2h/rPsdA
zxq7Pxyl0DJ7OVe37G9MGNfivvkkgeuJAueVLrUo2UocL8uSQPRNC08PWR0u
SAvKDPGQwG7WIemvjS9w23TlWi9tCSSruPkhRD3ON70xPCYchF07HXS/Br3E
1pvuFTdTBuFZuktE1IoGDEtat1PLZxBaQv21dIYa8LbOq+4Ki0EInyX3Mq8W
YsnmT9/jFQMgKPxwdXFiIwoMo6YFNQNwTmdXSq5fEyqc3assrw2AfLL54YZV
zajdmMj1YA3ADx/4XZXyZlxq4mVfupXeZzAeVNW8wmUF42Fs8wGYcKhOKLrY
glJjPXtySgzy0NSkYuZrbIhw0zR9EEP40dZNNqataHYynrKrEwPhd1p+vb8V
p/XT97HuiGHTwiteFyraUPfwmIcdIYYs4ew//WPbsfy2Z87+02K4qlXkzHR9
g/lLl+8IPygGk3n8+CqjDkz68aHIgCkGpqHxE8N3HXiHV2bNWS8GRuSkaNSv
E9VuQ72FFmJ4YH5GqqrpxO0xd2eCDcRQUc8b7bV6i+zzkc6sKQomO1IP6qe9
xaWmM3NshinIeWdVt0b1Ft8v2iT40k9BCq9+835vEbp6qphdHRRU/2O5sOWJ
CCMPZAYkvqTAXctOYbuExPk1iVxRFQWs9AiXQTMS90Zu0efRTltcknbLgsTW
jaqsk7Q711vZ6luTqPPUoEuP9qHQeWzZahJtRqgZeExBTMOYiudI4hbfrKHy
SgpKLtWOmOwl8Z1+9mR+OQUt86a2tvuSyMkZt2fTVhY4Zyb7kRjvKIpype1Q
Ve4wGkBicW2feW8ZBbWy4ujOYBLXbKdMTGl370wdS2OTqP4p5nL2X/R/soG5
LZLEmHFuJou25RFG9sQpEn36PwmcaIfHnnM6Hk2i5y47Rvd9Cr7wT/zOjCPx
2OPx3aa0jR35zd/jSRSs4NZ+LqXAqV6y6BGHRGm6GzyjndAT9HRZEoklpOrY
AdpFR2/p9V4gcf8Z2Rxb2i9HugIzkkkMzdMIp+5RIOUYCzxS6L/3dFLbaOvN
8pn4N5VEof5K/wLa/wH0WHWp
       "]]},
     Annotation[#, "Charting`Private`Tag$12966853#1"]& ]}, {}, {}},
  AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
  Axes->{True, True},
  AxesLabel->{None, None},
  AxesOrigin->{0, 0},
  DisplayFunction->Identity,
  Epilog->{
    PointBox[{1, 0.9815342232285977}], 
    PointBox[{10, 0.9837662804827835}], 
    PointBox[{50, 0.9789464510688851}], 
    PointBox[{100, 0.9763015390525981}], 
    PointBox[{500, 0.9452291302847118}], 
    PointBox[{1000, 0.9125586085549242}], 
    PointBox[{1500, 0.880238284394585}], 
    PointBox[{2000, 0.8555202736347571}], 
    PointBox[{3000, 0.800961364344768}], 
    PointBox[{5000, 0.7179087781596724}], 
    PointBox[{10000, 0.6024777955449644}], 
    PointBox[{20000, 0.520055771149573}]},
  Frame->{{False, False}, {False, False}},
  FrameLabel->{{None, None}, {None, None}},
  FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
  GridLines->{None, None},
  GridLinesStyle->Directive[
    GrayLevel[0.5, 0.4]],
  ImagePadding->All,
  Method->{
   "DefaultBoundaryStyle" -> Automatic, "DefaultMeshStyle" -> 
    AbsolutePointSize[6], "ScalingFunctions" -> None, 
    "CoordinatesToolOptions" -> {"DisplayFunction" -> ({
        (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
         Part[#, 1]], 
        (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
         Part[#, 2]]}& ), "CopiedValueFunction" -> ({
        (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
         Part[#, 1]], 
        (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
         Part[#, 2]]}& )}},
  PlotRange->{All, {0, 1}},
  PlotRangeClipping->True,
  PlotRangePadding->{{
     Scaled[0.02], 
     Scaled[0.02]}, {0, 0}},
  Ticks->{Automatic, Automatic}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"Histogram", "[", 
    RowBox[{"#", ",", 
     RowBox[{"{", 
      RowBox[{"0", ",", "1", ",", "0.01"}], "}"}], ",", 
     RowBox[{"ImageSize", "\[Rule]", "400"}]}], "]"}], "&"}], "/@", 
  RowBox[{"survivals", "[", 
   RowBox[{"[", 
    RowBox[{"All", ",", "1", ",", "All", ",", "1"}], "]"}], "]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.9500000000000001, 0}, {0.96, 9},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 28},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 21},
        RoundingRadius->0], RectangleBox[{0.99, 0}, {1., 28},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.9500000000000001, 0}, {0.96, 12},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 21},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 19},
        RoundingRadius->0], RectangleBox[{0.99, 0}, {1., 48},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 6},
        RoundingRadius->0], RectangleBox[{0.9500000000000001, 0}, {0.96, 16},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 15},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 26},
        RoundingRadius->0], RectangleBox[{0.99, 0}, {1., 37},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.93, 0}, {0.9400000000000001, 2},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 6},
        RoundingRadius->0], RectangleBox[{0.9500000000000001, 0}, {0.96, 17},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 9},
        RoundingRadius->0], RectangleBox[{0.97, 0}, {0.98, 9},
        RoundingRadius->0], RectangleBox[{0.98, 0}, {0.99, 21},
        RoundingRadius->0], RectangleBox[{0.99, 0}, {1., 36},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.84, 0}, {0.85, 1},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 1},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 1},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 6},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 4},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 7},
        RoundingRadius->0], RectangleBox[{0.93, 0}, {0.9400000000000001, 11},
        RoundingRadius->0], 
       RectangleBox[{0.9400000000000001, 0}, {0.9500000000000001, 15},
        RoundingRadius->0], RectangleBox[{0.9500000000000001, 0}, {0.96, 27},
        RoundingRadius->0], RectangleBox[{0.96, 0}, {0.97, 27},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 1},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 2},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 2},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 3},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 3},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 6},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 9},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 8},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 15},
        RoundingRadius->0], RectangleBox[{0.92, 0}, {0.93, 18},
        RoundingRadius->0], RectangleBox[{0.93, 0}, {0.9400000000000001, 33},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.76, 0}, {0.77, 1},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 1},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 2},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 2},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 1},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 1},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 5},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 7},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 9},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 8},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 11},
        RoundingRadius->0], RectangleBox[{0.89, 0}, {0.9, 20},
        RoundingRadius->0], RectangleBox[{0.9, 0}, {0.91, 31},
        RoundingRadius->0], RectangleBox[{0.91, 0}, {0.92, 1},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.72, 0}, {0.73, 1},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 1},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 1},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 2},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 2},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 2},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 6},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 4},
        RoundingRadius->0], RectangleBox[{0.84, 0}, {0.85, 10},
        RoundingRadius->0], RectangleBox[{0.85, 0}, {0.86, 12},
        RoundingRadius->0], RectangleBox[{0.86, 0}, {0.87, 21},
        RoundingRadius->0], RectangleBox[{0.87, 0}, {0.88, 28},
        RoundingRadius->0], RectangleBox[{0.88, 0}, {0.89, 10},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 1},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 2},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 1},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 3},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 2},
        RoundingRadius->0], RectangleBox[{0.76, 0}, {0.77, 2},
        RoundingRadius->0], RectangleBox[{0.77, 0}, {0.78, 8},
        RoundingRadius->0], RectangleBox[{0.78, 0}, {0.79, 12},
        RoundingRadius->0], RectangleBox[{0.79, 0}, {0.8, 10},
        RoundingRadius->0], RectangleBox[{0.8, 0}, {0.81, 10},
        RoundingRadius->0], RectangleBox[{0.81, 0}, {0.8200000000000001, 19},
        RoundingRadius->0], 
       RectangleBox[{0.8200000000000001, 0}, {0.8300000000000001, 20},
        RoundingRadius->0], RectangleBox[{0.8300000000000001, 0}, {0.84, 10},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.61, 0}, {0.62, 1},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 5},
        RoundingRadius->0], RectangleBox[{0.64, 0}, {0.65, 1},
        RoundingRadius->0], RectangleBox[{0.65, 0}, {0.66, 1},
        RoundingRadius->0], RectangleBox[{0.66, 0}, {0.67, 6},
        RoundingRadius->0], RectangleBox[{0.67, 0}, {0.68, 1},
        RoundingRadius->0], RectangleBox[{0.68, 0}, {0.6900000000000001, 2},
        RoundingRadius->0], 
       RectangleBox[{0.6900000000000001, 0}, {0.7000000000000001, 3},
        RoundingRadius->0], RectangleBox[{0.7000000000000001, 0}, {0.71, 6},
        RoundingRadius->0], RectangleBox[{0.71, 0}, {0.72, 14},
        RoundingRadius->0], RectangleBox[{0.72, 0}, {0.73, 12},
        RoundingRadius->0], RectangleBox[{0.73, 0}, {0.74, 19},
        RoundingRadius->0], RectangleBox[{0.74, 0}, {0.75, 16},
        RoundingRadius->0], RectangleBox[{0.75, 0}, {0.76, 13},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, \
{}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.52, 0}, {0.53, 2},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 2},
        RoundingRadius->0], RectangleBox[{0.54, 0}, {0.55, 1},
        RoundingRadius->0], RectangleBox[{0.55, 0}, {0.56, 6},
        RoundingRadius->0], RectangleBox[{0.56, 0}, {0.5700000000000001, 5},
        RoundingRadius->0], RectangleBox[{0.5700000000000001, 0}, {0.58, 3},
        RoundingRadius->0], RectangleBox[{0.58, 0}, {0.59, 5},
        RoundingRadius->0], RectangleBox[{0.59, 0}, {0.6, 10},
        RoundingRadius->0], RectangleBox[{0.6, 0}, {0.61, 17},
        RoundingRadius->0], RectangleBox[{0.61, 0}, {0.62, 19},
        RoundingRadius->0], RectangleBox[{0.62, 0}, {0.63, 23},
        RoundingRadius->0], RectangleBox[{0.63, 0}, {0.64, 7},
        RoundingRadius->
         0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}], ",", 
   GraphicsBox[{
     {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
      Opacity[0.]], {}, 
      {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[
       Opacity[0.]], RectangleBox[{0.47000000000000003, 0}, {0.48, 1},
        RoundingRadius->0], RectangleBox[{0.48, 0}, {0.49, 1},
        RoundingRadius->0], RectangleBox[{0.49, 0}, {0.5, 8},
        RoundingRadius->0], RectangleBox[{0.5, 0}, {0.51, 12},
        RoundingRadius->0], RectangleBox[{0.51, 0}, {0.52, 20},
        RoundingRadius->0], RectangleBox[{0.52, 0}, {0.53, 32},
        RoundingRadius->0], RectangleBox[{0.53, 0}, {0.54, 26},
        RoundingRadius->0]}, {}, {}}, {{}, {}, {}, {}, {}, {}, {}}},
    AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
    Axes->{True, True},
    AxesLabel->{None, None},
    AxesOrigin->{-0.02, 0},
    FrameLabel->{{None, None}, {None, None}},
    FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
    GridLines->{None, None},
    GridLinesStyle->Directive[
      GrayLevel[0.5, 0.4]],
    ImageSize->400,
    PlotRange->{{0., 1.}, {All, All}},
    PlotRangePadding->{{
       Scaled[0.02], 
       Scaled[0.02]}, {
       Scaled[0.02], 
       Scaled[0.05]}},
    Ticks->{Automatic, Automatic}]}], "}"}]], "Output"]
}, Open  ]]
}, Closed]],

Cell[CellGroupData[{

Cell["For LRB", "Section"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"Dij", "[", 
   RowBox[{"i_", ",", "j_", ",", "d1_", ",", "d2_"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"Ii", ",", 
      RowBox[{"d", "=", 
       RowBox[{"d1", "+", "d2"}]}]}], "}"}], ",", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"If", "[", 
      RowBox[{
       RowBox[{"i", "\[Equal]", "1"}], ",", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"Ii", "=", 
         RowBox[{"DiagonalMatrix", "[", 
          RowBox[{"Join", "[", 
           RowBox[{
            RowBox[{"ConstantArray", "[", 
             RowBox[{
              RowBox[{"1", "/", "d1"}], ",", "d1"}], "]"}], ",", 
            RowBox[{"ConstantArray", "[", 
             RowBox[{"0", ",", "d2"}], "]"}]}], "]"}], "]"}]}], ";"}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"Ii", "=", 
         RowBox[{"DiagonalMatrix", "[", 
          RowBox[{"Join", "[", 
           RowBox[{
            RowBox[{"ConstantArray", "[", 
             RowBox[{"0", ",", "d1"}], "]"}], ",", 
            RowBox[{"ConstantArray", "[", 
             RowBox[{
              RowBox[{"1", "/", "d2"}], ",", "d2"}], "]"}]}], "]"}], "]"}]}], 
        ";"}]}], "\[IndentingNewLine]", "]"}], ";", "\[IndentingNewLine]", 
     RowBox[{"If", "[", 
      RowBox[{
       RowBox[{"j", "\[Equal]", "1"}], ",", "\[IndentingNewLine]", 
       RowBox[{"Super", "@", 
        RowBox[{"FunctionChannel", "[", 
         RowBox[{
          RowBox[{
           RowBox[{
            RowBox[{"Total", "[", 
             RowBox[{
              RowBox[{"Diagonal", "[", "#", "]"}], "[", 
              RowBox[{"[", 
               RowBox[{";;", "d1"}], "]"}], "]"}], "]"}], "Ii"}], "&"}], ",", 
          
          RowBox[{"InputDim", "->", "d"}]}], "]"}]}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{"Super", "@", 
        RowBox[{"FunctionChannel", "[", 
         RowBox[{
          RowBox[{
           RowBox[{
            RowBox[{"Total", "[", 
             RowBox[{
              RowBox[{"Diagonal", "[", "#", "]"}], "[", 
              RowBox[{"[", 
               RowBox[{
                RowBox[{"d1", "+", "1"}], ";;"}], "]"}], "]"}], "]"}], "Ii"}],
            "&"}], ",", 
          RowBox[{"InputDim", "->", "d"}]}], "]"}]}]}], "\[IndentingNewLine]",
       "]"}]}]}], "\[IndentingNewLine]", "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"L", "[", 
   RowBox[{"\[Rho]_", ",", "d1_"}], "]"}], ":=", 
  RowBox[{"1", "-", 
   RowBox[{"Total", "[", 
    RowBox[{
     RowBox[{"Diagonal", "[", "\[Rho]", "]"}], "[", 
     RowBox[{"[", 
      RowBox[{";;", "d1"}], "]"}], "]"}], "]"}]}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"L1", "[", 
   RowBox[{"E_", ",", "d1_"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"I1", ",", 
      RowBox[{"d2", "=", 
       RowBox[{
        RowBox[{"InputDim", "@", "E"}], "-", "d1"}]}]}], "}"}], ",", 
    "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"I1", "=", 
      RowBox[{"DiagonalMatrix", "[", 
       RowBox[{"Join", "[", 
        RowBox[{
         RowBox[{"ConstantArray", "[", 
          RowBox[{"1", ",", "d1"}], "]"}], ",", 
         RowBox[{"ConstantArray", "[", 
          RowBox[{"0", ",", "d2"}], "]"}]}], "]"}], "]"}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"L", "[", 
      RowBox[{
       RowBox[{"E", "[", 
        RowBox[{"I1", "/", "d1"}], "]"}], ",", "d1"}], "]"}]}]}], 
   "\[IndentingNewLine]", "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"L2", "[", 
   RowBox[{"E_", ",", "d1_"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"I2", ",", 
      RowBox[{"d2", "=", 
       RowBox[{
        RowBox[{"InputDim", "@", "E"}], "-", "d1"}]}]}], "}"}], ",", 
    "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"I2", "=", 
      RowBox[{"DiagonalMatrix", "[", 
       RowBox[{"Join", "[", 
        RowBox[{
         RowBox[{"ConstantArray", "[", 
          RowBox[{"0", ",", "d1"}], "]"}], ",", 
         RowBox[{"ConstantArray", "[", 
          RowBox[{"1", ",", "d2"}], "]"}]}], "]"}], "]"}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"1", "-", 
      RowBox[{"L", "[", 
       RowBox[{
        RowBox[{"E", "[", 
         RowBox[{"I2", "/", "d2"}], "]"}], ",", "d1"}], "]"}]}]}]}], 
   "\[IndentingNewLine]", "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"F1", "[", 
   RowBox[{"E_", ",", "d1_"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"d2", "=", 
       RowBox[{
        RowBox[{"InputDim", "@", "E"}], "-", "d1"}]}], ",", "I1", ",", 
      "Fpro"}], "}"}], ",", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"I1", "=", 
      RowBox[{"Join", "[", 
       RowBox[{
        RowBox[{"ConstantArray", "[", 
         RowBox[{"1", ",", "d1"}], "]"}], ",", 
        RowBox[{"ConstantArray", "[", 
         RowBox[{"0", ",", "d2"}], "]"}]}], "]"}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"Fpro", "=", 
      RowBox[{
       RowBox[{"Total", "[", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{"I1", "\[CircleTimes]", "I1"}], ")"}], "*", 
         RowBox[{"Diagonal", "[", 
          RowBox[{"First", "@", 
           RowBox[{"Super", "@", "E"}]}], "]"}]}], "]"}], "/", 
       RowBox[{"d1", "^", "2"}]}]}], ";", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{
        RowBox[{"d1", " ", "Fpro"}], "+", "1", "-", 
        RowBox[{"L1", "[", 
         RowBox[{"E", ",", "d1"}], "]"}]}], ")"}], "/", 
      RowBox[{"(", 
       RowBox[{"d1", "+", "1"}], ")"}]}]}]}], "\[IndentingNewLine]", 
   "]"}]}]}], "Input"],

Cell["\<\
Here is the depolarizing leakage extension (DLE) for a map E acting only on \
space 1.\
\>", "Text"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"AddIdentity", "[", 
   RowBox[{"K_List", ",", "dim_"}], "]"}], ":=", 
  RowBox[{"BlockMatrix", "[", 
   RowBox[{"K", ",", 
    RowBox[{"0", "*", 
     RowBox[{"IdentityMatrix", "[", "dim", "]"}]}]}], 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AddIdentity", "[", 
   RowBox[{"E_", ",", "dim_"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"all", "=", 
       RowBox[{"First", "@", 
        RowBox[{"Kraus", "[", "E", "]"}]}]}], ",", "Ks", ",", "Ls"}], "}"}], 
    ",", "\[IndentingNewLine]", 
    RowBox[{"Super", "@", 
     RowBox[{"If", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"Length", "@", 
         RowBox[{"Dimensions", "@", "all"}]}], "\[Equal]", "3"}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{"Kraus", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"AddIdentity", "[", 
           RowBox[{"#", ",", "dim"}], "]"}], "&"}], "/@", "all"}], "]"}], ",",
        "\[IndentingNewLine]", 
       RowBox[{"Kraus", "[", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           RowBox[{
            RowBox[{"AddIdentity", "[", 
             RowBox[{"#", ",", "dim"}], "]"}], "&"}], "/@", 
           RowBox[{"all", "[", 
            RowBox[{"[", "1", "]"}], "]"}]}], ",", 
          RowBox[{
           RowBox[{
            RowBox[{"AddIdentity", "[", 
             RowBox[{"#", ",", "dim"}], "]"}], "&"}], "/@", 
           RowBox[{"all", "[", 
            RowBox[{"[", "2", "]"}], "]"}]}]}], "}"}], "]"}]}], 
      "\[IndentingNewLine]", "]"}]}]}], "\[IndentingNewLine]", 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"DLE", "[", 
   RowBox[{"E_", ",", "L1_", ",", "L2_", ",", "d2_"}], "]"}], ":=", 
  RowBox[{"With", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"d1", "=", 
      RowBox[{"InputDim", "@", "E"}]}], "}"}], ",", 
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1", "-", "L1"}], ")"}], 
      RowBox[{"AddIdentity", "[", 
       RowBox[{"E", ",", "d2"}], "]"}]}], "+", 
     RowBox[{"L1", " ", 
      RowBox[{"Dij", "[", 
       RowBox[{"2", ",", "1", ",", "d1", ",", "d2"}], "]"}]}], "+", 
     RowBox[{"L2", " ", 
      RowBox[{"Dij", "[", 
       RowBox[{"1", ",", "2", ",", "d1", ",", "d2"}], "]"}]}], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1", "-", "L2"}], ")"}], 
      RowBox[{"Dij", "[", 
       RowBox[{"2", ",", "2", ",", "d1", ",", "d2"}], "]"}]}]}]}], 
   "]"}]}]}], "Input"],

Cell["Now set up a gate set with gate independent noise", "Text"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"dephasingGate", "[", "p_", "]"}], ":=", 
  RowBox[{"Kraus", "[", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{
      RowBox[{"Sqrt", "[", 
       RowBox[{"1", "-", "p"}], "]"}], 
      RowBox[{"TP", "[", "I", "]"}]}], ",", 
     RowBox[{
      RowBox[{"Sqrt", "[", "p", "]"}], 
      RowBox[{"TP", "[", "Z", "]"}]}]}], "}"}], 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{"unitaryGate", ":=", 
  RowBox[{"Unitary", "[", 
   RowBox[{"MatrixExp", "[", 
    RowBox[{
     RowBox[{"-", "I"}], " ", 
     RowBox[{"(", 
      RowBox[{"0.1", 
       RowBox[{"\[Pi]", "/", "180"}]}], ")"}], 
     RowBox[{
      RowBox[{"TP", "[", "Z", "]"}], "/", "2"}]}], "]"}], 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateNoise", "=", 
   RowBox[{"DLE", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"dephasingGate", "[", "0.003", "]"}], ".", "unitaryGate"}], ",",
      "0.001", ",", "0.0015", ",", "1"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Re", "@", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"F1", "[", 
      RowBox[{"gateNoise", ",", "2"}], "]"}], ",", 
     RowBox[{"L1", "[", 
      RowBox[{"gateNoise", ",", "2"}], "]"}], ",", 
     RowBox[{"L2", "[", 
      RowBox[{"gateNoise", ",", "2"}], "]"}]}], "}"}]}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateUnitaries", "=", 
   RowBox[{
    RowBox[{
     RowBox[{"BlockMatrix", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"GateUnitary", "[", "$gaussianQubitGateSet", "]"}], "[", "#", 
        "]"}], ",", 
       RowBox[{"{", 
        RowBox[{"{", "1", "}"}], "}"}]}], "]"}], "&"}], "/@", 
    RowBox[{"Range", "[", "12", "]"}]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"gateChannels", "=", 
    RowBox[{
     RowBox[{
      RowBox[{"Super", "[", 
       RowBox[{
        RowBox[{"Unitary", "[", "#", "]"}], ".", "gateNoise"}], "]"}], "&"}], 
     "/@", "gateUnitaries"}]}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gs", "=", 
   RowBox[{"GateNoiseGateSet", "[", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
     "GateSetName", "\[Rule]", "\"\<Z+Dephasing with L1 and L2 leakage\>\""}],
      ",", "\[IndentingNewLine]", 
     RowBox[{"Size", "\[Rule]", "12"}], ",", "\[IndentingNewLine]", 
     RowBox[{"Dimension", "\[Rule]", "2"}], ",", "\[IndentingNewLine]", 
     RowBox[{"GateProduct", "\[Rule]", 
      RowBox[{"GateProduct", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateInverse", "\[Rule]", 
      RowBox[{"GateInverse", "[", "$gaussianQubitGateSet", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateUnitary", "\[Rule]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"gateUnitaries", "[", 
         RowBox[{"[", 
          RowBox[{"Last", "[", 
           RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}], ",",
      "\[IndentingNewLine]", 
     RowBox[{"GateNoise", "\[Rule]", 
      RowBox[{"IndependentNoise", "[", "gateNoise", "]"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"GateChannel", "\[Rule]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"gateChannels", "[", 
         RowBox[{"[", 
          RowBox[{"Last", "[", 
           RowBox[{"{", "##", "}"}], "]"}], "]"}], "]"}], "&"}], ")"}]}], ",",
      "\[IndentingNewLine]", 
     RowBox[{"GateNoiseMemoryDepth", "\[Rule]", "1"}]}], 
    "\[IndentingNewLine]", "]"}]}], ";"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0.9970014958552521`", ",", "0.0010000000000000009`", ",", 
   "0.0014999999999999458`"}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"LRBProtocol", "[", 
   RowBox[{
   "seqLengths_", ",", "numSeqs_", ",", "shotsPerSeq_", ",", "\[Rho]list_", 
    ",", "Mlist_"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"seqGen", ",", "gateSim"}], "}"}], ",", "\n", "\t", 
    RowBox[{"(*", " ", 
     RowBox[{
     "Use", " ", "two", " ", "closures", " ", "to", " ", "make", " ", "sure", 
      " ", "every", " ", "shot", " ", "of", " ", "the", " ", "same", " ", 
      RowBox[{"(", 
       RowBox[{"seqLength", ",", "seqNum"}], ")"}], " ", "is", " ", 
      "identical"}], " ", "*)"}], "\n", "\t", 
    RowBox[{
     RowBox[{
      RowBox[{"seqGen", "[", 
       RowBox[{"gs_", ",", 
        RowBox[{"{", 
         RowBox[{"seqLength_", ",", "_", ",", "seqNum_"}], "}"}]}], "]"}], 
      " ", ":=", 
      RowBox[{"RBDraw", "[", 
       RowBox[{"gs", ",", " ", "seqLength"}], "]"}]}], ";", "\n", "\t", 
     RowBox[{
      RowBox[{"gateSim", "[", 
       RowBox[{"gs_GateNoiseGateSet", ",", 
        RowBox[{"{", "idxs__", "}"}]}], "]"}], " ", ":=", " ", 
      RowBox[{"With", "[", 
       RowBox[{
        RowBox[{"{", "\[IndentingNewLine]", 
         RowBox[{"op", "=", 
          RowBox[{"Fold", "[", "\[IndentingNewLine]", 
           RowBox[{
            RowBox[{
             RowBox[{"(", 
              RowBox[{
               RowBox[{"First", "[", 
                RowBox[{"Super", "[", 
                 RowBox[{
                  RowBox[{"GateChannel", "[", "gs", "]"}], "[", "#2", "]"}], 
                 "]"}], "]"}], ".", "#1"}], ")"}], "&"}], ",", 
            "\[IndentingNewLine]", 
            RowBox[{"Join", "[", 
             RowBox[{
              RowBox[{"Sequence", "@@", 
               RowBox[{"(", 
                RowBox[{"Vec", "/@", "\[Rho]list"}], ")"}]}], ",", "2"}], 
             "]"}], ",", "\[IndentingNewLine]", " ", 
            RowBox[{"{", "idxs", "}"}]}], "\[IndentingNewLine]", "]"}]}], 
         "}"}], ",", "\[IndentingNewLine]", 
        RowBox[{"Re", "[", 
         RowBox[{
          RowBox[{"(", 
           RowBox[{
            RowBox[{
             RowBox[{"Flatten", "[", 
              RowBox[{"Vec", "[", "#", "]"}], "]"}], "&"}], "/@", "Mlist"}], 
           ")"}], ".", "op"}], "]"}]}], "\[IndentingNewLine]", "]"}]}], ";", 
     "\n", "\t", 
     RowBox[{"Protocol", "[", "\n", "\t\t", 
      RowBox[{
       RowBox[{
       "ProtocolName", " ", "->", " ", "\"\<Randomized Benchmarking\>\""}], 
       ",", "\n", "\t\t", 
       RowBox[{"SequenceLengths", " ", "->", " ", "seqLengths"}], ",", "\n", 
       "\t\t", 
       RowBox[{"ExperimentTypes", " ", "->", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"{", "0", "}"}], "&"}], ")"}]}], ",", "\n", "\t\t", 
       RowBox[{"NumSequenceDraws", " ", "->", " ", 
        RowBox[{"(", 
         RowBox[{"numSeqs", "&"}], ")"}]}], ",", "\n", "\t\t", 
       RowBox[{"NumRepetitions", " ", "->", " ", 
        RowBox[{"(", 
         RowBox[{"shotsPerSeq", "&"}], ")"}]}], ",", "\n", "\t\t", 
       RowBox[{"SequenceGenerator", " ", "->", " ", "seqGen"}], ",", "\n", 
       "\t\t", 
       RowBox[{"SimulationOptions", " ", "->", " ", 
        RowBox[{"NotImplemented", "[", "]"}]}], ",", "\n", "\t\t", 
       RowBox[{"SimulationParser", " ", "->", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           RowBox[{"Observables", "[", "#", "]"}], "[", 
           RowBox[{"[", 
            RowBox[{"All", ",", 
             RowBox[{"-", "1"}]}], "]"}], "]"}], "&"}], ")"}]}], ",", "\n", 
       "\t\t", 
       RowBox[{"GateSimulator", " ", "->", " ", "gateSim"}], ",", "\n", 
       "\t\t", 
       RowBox[{"ParallelOptions", " ", "->", " ", 
        RowBox[{"{", 
         RowBox[{"Method", "->", "\"\<CoarsestGrained\>\""}], "}"}]}], ",", 
       "\n", "\t\t", 
       RowBox[{"TotalGates", " ", "->", " ", 
        RowBox[{"(", 
         RowBox[{"shotsPerSeq", " ", "*", " ", "numSeqs", " ", "*", " ", 
          RowBox[{"Total", "[", 
           RowBox[{"seqLengths", "+", "1"}], "]"}]}], ")"}]}]}], "\n", "\t", 
      "]"}]}]}], "\n", "]"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"mlist", "=", 
   RowBox[{"{", 
    RowBox[{
    "1", ",", "10", ",", "50", ",", "100", ",", "200", ",", "300", ",", "500",
      ",", "1000", ",", "1500", ",", "2000", ",", "3000", ",", "4000"}], 
    "}"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"protocol", "=", 
   RowBox[{"LRBProtocol", "[", 
    RowBox[{"mlist", ",", "100", ",", "1", ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"Projector", "[", 
        RowBox[{"{", 
         RowBox[{"0.9999", ",", "0", ",", "0"}], "}"}], "]"}], ",", 
       RowBox[{"Projector", "[", 
        RowBox[{"{", 
         RowBox[{"0", ",", "0.9995", ",", "0"}], "}"}], "]"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"Projector", "[", 
        RowBox[{"{", 
         RowBox[{"0.99999", ",", "0", ",", "0"}], "}"}], "]"}], ",", 
       RowBox[{"Projector", "[", 
        RowBox[{"{", 
         RowBox[{"0", ",", "0.999", ",", "0"}], "}"}], "]"}]}], "}"}]}], 
    "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"survivals", "=", 
   RowBox[{"SimulateProtocol", "[", 
    RowBox[{"gs", ",", "protocol", ",", 
     RowBox[{
     "SimulationExportName", "\[Rule]", 
      "\"\<../data/lrb-dephasing_z_andL1L2\>\""}]}], "]"}]}], ";"}]}], "Input"],

Cell[BoxData["\<\"Saved ../data/lrb-dephasing_z_andL1L2.h5\"\>"], "Print"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"survivals", "//", "Dimensions"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"9", ",", "1", ",", "500", ",", "1", ",", "2", ",", "2"}], 
  "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"ListPlot", "[", 
  RowBox[{"{", "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     RowBox[{"{", 
      RowBox[{"mlist", ",", 
       RowBox[{"Mean", "/@", 
        RowBox[{"survivals", "[", 
         RowBox[{"[", 
          RowBox[{"All", ",", "1", ",", "All", ",", "1", ",", "1", ",", "1"}],
           "]"}], "]"}]}]}], "}"}], "\[Transpose]"}], ",", 
    "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"mlist", ",", 
       RowBox[{"Mean", "/@", 
        RowBox[{"survivals", "[", 
         RowBox[{"[", 
          RowBox[{"All", ",", "1", ",", "All", ",", "1", ",", "1", ",", "2"}],
           "]"}], "]"}]}]}], "}"}], "\[Transpose]"}]}], "\[IndentingNewLine]",
    "}"}], "]"}]], "Input"],

Cell[BoxData[
 GraphicsBox[{{}, {{}, 
    {RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.012833333333333334`],
      AbsoluteThickness[1.6], 
     PointBox[{{1., 0.9958994295604392}, {50., 0.864391556219469}, {100., 
      0.7583771485350492}, {200., 0.6045515668647682}, {300., 
      0.5051086739482625}, {400., 0.44089964840696905`}, {500., 
      0.3979507716862442}, {1000., 0.3196544918024357}, {1500., 
      0.3049431233952251}, {2000., 0.3013582256037278}, {3000., 
      0.30010945703813885`}, {4000., 0.3000089455745027}}]}, 
    {RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.012833333333333334`],
      AbsoluteThickness[1.6], 
     PointBox[{{1., 0.002103070439555221}, {50., 0.0876696021887294}, {100., 
      0.1522671102612482}, {200., 0.23730270762921912`}, {300., 
      0.28318867401389586`}, {400., 0.30570057738504486`}, {500., 
      0.31618586610957056`}, {1000., 0.3129950623011603}, {1500., 
      0.30439650149429925`}, {2000., 0.301313437010732}, {3000., 
      0.3001091609896506}, {4000., 
      0.30000894360129127`}}]}, {}}, {}, {}, {{}, {}}},
  AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
  Axes->{True, True},
  AxesLabel->{None, None},
  AxesOrigin->{0, 0},
  DisplayFunction->Identity,
  Frame->{{False, False}, {False, False}},
  FrameLabel->{{None, None}, {None, None}},
  FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
  GridLines->{None, None},
  GridLinesStyle->Directive[
    GrayLevel[0.5, 0.4]],
  ImagePadding->All,
  Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
        (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
         Part[#, 1]], 
        (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
         Part[#, 2]]}& ), "CopiedValueFunction" -> ({
        (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
         Part[#, 1]], 
        (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
         Part[#, 2]]}& )}},
  PlotRange->{{0, 4000.}, {0, 0.9958994295604392}},
  PlotRangeClipping->True,
  PlotRangePadding->{{
     Scaled[0.02], 
     Scaled[0.02]}, {
     Scaled[0.02], 
     Scaled[0.05]}},
  Ticks->{Automatic, Automatic}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{
    RowBox[{"MatrixPower", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"GateChannel", "[", "gs", "]"}], "[", "1", "]"}], ",", 
      "50000"}], "]"}], "[", 
    RowBox[{"DiagonalMatrix", "[", 
     RowBox[{"{", 
      RowBox[{"0", ",", "0", ",", "1"}], "}"}], "]"}], "]"}], "//", "Chop"}], 
  "//", "MatrixForm"}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"0.30000000000075766`", "0", "0"},
     {"0", "0.30000000000075766`", "0"},
     {"0", "0", "0.400000000001014`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{"GateChannel", "[", "gs", "]"}], "[", "6", "]"}], "//", 
  "TracePreservingQ"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"dephasingGate", "[", "0.001", "]"}], "//", 
  "TracePreservingQ"}]}], "Input"],

Cell[BoxData["True"], "Output"],

Cell[BoxData["True"], "Output"]
}, Open  ]]
}, Closed]],

Cell[CellGroupData[{

Cell["Nonmarkovian", "Section"],

Cell[BoxData[
 RowBox[{
  RowBox[{"nmm", "[", 
   RowBox[{"J_", ",", "T2_", ",", "T2env_"}], "]"}], ":=", 
  RowBox[{"NoiseModel", "[", "\[IndentingNewLine]", 
   RowBox[{"Generator", "\[Rule]", 
    RowBox[{"LindbladForm", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"2", "\[Pi]", " ", "J", " ", 
        RowBox[{
         RowBox[{"TP", "[", "ZZ", "]"}], "/", "4"}]}], "+", 
       RowBox[{"TP", "[", "IX", "]"}]}], ",", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"Sqrt", "[", 
          RowBox[{"1", "/", "T2"}], "]"}], 
         RowBox[{"TP", "[", "ZI", "]"}]}], ",", 
        RowBox[{
         RowBox[{"Sqrt", "[", 
          RowBox[{"1", "/", "T2env"}], "]"}], 
         RowBox[{"TP", "[", "IZ", "]"}]}]}], "}"}]}], "]"}]}], 
   "\[IndentingNewLine]", "]"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"rawPulseX2", "=", 
   RowBox[{"GaussianTailsPulse", "[", 
    RowBox[{"0.001", ",", "0.01", ",", "0.005", ",", 
     RowBox[{"Area", "\[Rule]", 
      RowBox[{"1", "/", "4"}]}]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"rawPulseY2", "=", 
   RowBox[{"FromPulse", "@", 
    RowBox[{"PulsePhaseRotate", "[", 
     RowBox[{
      RowBox[{"ToPulse", "@", "rawPulseX2"}], ",", 
      RowBox[{"\[Pi]", "/", "2"}]}], "]"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"rawPulseXM2", "=", 
   RowBox[{"FromPulse", "@", 
    RowBox[{"PulsePhaseRotate", "[", 
     RowBox[{
      RowBox[{"ToPulse", "@", "rawPulseX2"}], ",", "\[Pi]"}], "]"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"rawPulseYM2", "=", 
   RowBox[{"FromPulse", "@", 
    RowBox[{"PulsePhaseRotate", "[", 
     RowBox[{
      RowBox[{"ToPulse", "@", "rawPulseX2"}], ",", 
      RowBox[{
       RowBox[{"-", "\[Pi]"}], "/", "2"}]}], "]"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"group", "=", 
   RowBox[{"List", "@@@", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"Reverse", "/@", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"W", "[", "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"X2", ",", "X2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"X2", ",", "Y2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"XM2", ",", "Y2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"XM2", ",", "YM2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"YM2", ",", "XM2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"X2", ",", "YM2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"Y2", ",", "X2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"Y2", ",", "XM2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"YM2", ",", "X2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"Y2", ",", "Y2"}], "]"}], ",", 
         RowBox[{"W", "[", 
          RowBox[{"X2", ",", "X2", ",", "Y2", ",", "Y2"}], "]"}]}], "}"}]}], "/.", 
      RowBox[{
       RowBox[{"W", "[", "]"}], "\[Rule]", 
       RowBox[{"W", "[", "Id", "]"}]}]}], ")"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gatePulses", "=", 
   RowBox[{
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"Join", "@@", 
        RowBox[{"Reverse", "[", "#", "]"}]}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"2", "\[Pi]", " ", 
          RowBox[{
           RowBox[{"TP", "[", "X", "]"}], "/", "2"}]}], ",", 
         RowBox[{"2", "\[Pi]", " ", 
          RowBox[{
           RowBox[{"TP", "[", "Y", "]"}], "/", "2"}]}]}], "}"}]}], "}"}], 
     "&"}], "/@", 
    RowBox[{"(", 
     RowBox[{"group", "/.", 
      RowBox[{"{", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"X2", "\[Rule]", "rawPulseX2"}], ",", "\[IndentingNewLine]", 
        RowBox[{"Y2", "\[Rule]", "rawPulseY2"}], ",", "\[IndentingNewLine]", 
        RowBox[{"XM2", "\[Rule]", "rawPulseXM2"}], ",", "\[IndentingNewLine]", 
        RowBox[{"YM2", "\[Rule]", "rawPulseYM2"}], ",", "\[IndentingNewLine]", 
        RowBox[{"Id", "\[Rule]", 
         RowBox[{"{", 
          RowBox[{"{", 
           RowBox[{"$MachineEpsilon", ",", "0", ",", "0"}], "}"}], "}"}]}]}], 
       "\[IndentingNewLine]", "}"}]}], ")"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"gateUnitaries", "=", 
   RowBox[{"Dot", "@@@", 
    RowBox[{"(", 
     RowBox[{"group", "/.", 
      RowBox[{"{", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"X2", "->", 
         RowBox[{"MatrixExp", "[", 
          RowBox[{
           RowBox[{"-", "I"}], " ", 
           RowBox[{"(", 
            RowBox[{"\[Pi]", "/", "2"}], ")"}], 
           RowBox[{
            RowBox[{"TP", "[", "XI", "]"}], "/", "2"}]}], "]"}]}], ",", 
        "\[IndentingNewLine]", 
        RowBox[{"Y2", "->", 
         RowBox[{"MatrixExp", "[", 
          RowBox[{
           RowBox[{"-", "I"}], " ", 
           RowBox[{"(", 
            RowBox[{"\[Pi]", "/", "2"}], ")"}], 
           RowBox[{
            RowBox[{"TP", "[", "YI", "]"}], "/", "2"}]}], "]"}]}], ",", 
        "\[IndentingNewLine]", 
        RowBox[{"XM2", "->", 
         RowBox[{"MatrixExp", "[", 
          RowBox[{"I", " ", 
           RowBox[{"(", 
            RowBox[{"\[Pi]", "/", "2"}], ")"}], 
           RowBox[{
            RowBox[{"TP", "[", "XI", "]"}], "/", "2"}]}], "]"}]}], ",", 
        "\[IndentingNewLine]", 
        RowBox[{"YM2", "->", 
         RowBox[{"MatrixExp", "[", 
          RowBox[{"I", " ", 
           RowBox[{"(", 
            RowBox[{"\[Pi]", "/", "2"}], ")"}], 
           RowBox[{
            RowBox[{"TP", "[", "YI", "]"}], "/", "2"}]}], "]"}]}], ",", 
        "\[IndentingNewLine]", 
        RowBox[{"Id", "\[Rule]", 
         RowBox[{"IdentityMatrix", "[", "4", "]"}]}]}], "\[IndentingNewLine]",
        "}"}]}], ")"}]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"mulTable", "=", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
      "1", ",", "2", ",", "3", ",", "4", ",", "5", ",", "6", ",", "7", ",", 
       "8", ",", "9", ",", "10", ",", "11", ",", "12"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "2", ",", "1", ",", "4", ",", "3", ",", "7", ",", "8", ",", "5", ",", 
       "6", ",", "10", ",", "9", ",", "12", ",", "11"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "3", ",", "5", ",", "6", ",", "9", ",", "10", ",", "1", ",", "8", ",", 
       "11", ",", "12", ",", "2", ",", "7", ",", "4"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "4", ",", "7", ",", "8", ",", "10", ",", "9", ",", "2", ",", "6", ",", 
       "12", ",", "11", ",", "1", ",", "5", ",", "3"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "5", ",", "3", ",", "9", ",", "6", ",", "8", ",", "11", ",", "10", ",", 
       "1", ",", "2", ",", "12", ",", "4", ",", "7"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "6", ",", "10", ",", "1", ",", "12", ",", "2", ",", "3", ",", "11", ",",
        "7", ",", "4", ",", "5", ",", "8", ",", "9"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "7", ",", "4", ",", "10", ",", "8", ",", "6", ",", "12", ",", "9", ",", 
       "2", ",", "1", ",", "11", ",", "3", ",", "5"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "8", ",", "9", ",", "2", ",", "11", ",", "1", ",", "4", ",", "12", ",", 
       "5", ",", "3", ",", "7", ",", "6", ",", "10"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "9", ",", "8", ",", "11", ",", "2", ",", "12", ",", "5", ",", "1", ",", 
       "4", ",", "7", ",", "3", ",", "10", ",", "6"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "10", ",", "6", ",", "12", ",", "1", ",", "11", ",", "7", ",", "2", ",",
        "3", ",", "5", ",", "4", ",", "9", ",", "8"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "11", ",", "12", ",", "5", ",", "7", ",", "3", ",", "9", ",", "4", ",", 
       "10", ",", "6", ",", "8", ",", "1", ",", "2"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "12", ",", "11", ",", "7", ",", "5", ",", "4", ",", "10", ",", "3", ",",
        "9", ",", "8", ",", "6", ",", "2", ",", "1"}], "}"}]}], "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"invTable", "=", 
   RowBox[{"{", 
    RowBox[{
    "1", ",", "2", ",", "6", ",", "10", ",", "8", ",", "3", ",", "9", ",", 
     "5", ",", "7", ",", "4", ",", "11", ",", "12"}], "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"gs", "=", 
  RowBox[{"GateSet", "[", "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"GateSetName", "\[Rule]", "\"\<coupled to second qubit\>\""}], 
    ",", "\[IndentingNewLine]", 
    RowBox[{"Size", "\[Rule]", "12"}], ",", "\[IndentingNewLine]", 
    RowBox[{"Dimension", "\[Rule]", "4"}], ",", "\[IndentingNewLine]", 
    RowBox[{"GateProduct", "\[Rule]", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"mulTable", "[", 
        RowBox[{"[", 
         RowBox[{"#1", ",", "#2"}], "]"}], "]"}], "&"}], ")"}]}], ",", 
    "\[IndentingNewLine]", 
    RowBox[{"GateInverse", "\[Rule]", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"invTable", "[", 
        RowBox[{"[", "#", "]"}], "]"}], "&"}], ")"}]}], ",", 
    "\[IndentingNewLine]", 
    RowBox[{"GateUnitary", "\[Rule]", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"gateUnitaries", "[", 
        RowBox[{"[", "#", "]"}], "]"}], "&"}], ")"}]}], ",", 
    "\[IndentingNewLine]", 
    RowBox[{"GatePulse", "\[Rule]", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"gatePulses", "[", 
        RowBox[{"[", "#", "]"}], "]"}], "&"}], ")"}]}]}], 
   "\[IndentingNewLine]", "]"}]}]}], "Input"],

Cell[BoxData[
 RowBox[{"\<\"GateSet\"\>", "[", 
  RowBox[{
   RowBox[{
   "GateSetName", "\[Rule]", "\<\"\\\"coupled to second qubit\\\"\"\>"}], ",", 
   RowBox[{"Dimension", "\[Rule]", "4"}], ",", 
   RowBox[{"Size", "\[Rule]", "12"}], ",", "\<\"...\"\>"}], "]"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"mlist", "=", 
   RowBox[{"Round", "[", 
    RowBox[{"10", "^", 
     RowBox[{"Table", "[", 
      RowBox[{"m", ",", 
       RowBox[{"{", 
        RowBox[{"m", ",", "0", ",", "3", ",", "0.3"}], "}"}]}], "]"}]}], 
    "]"}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData["mlist"], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "1", ",", "2", ",", "4", ",", "8", ",", "16", ",", "32", ",", "63", ",", 
   "126", ",", "251", ",", "501", ",", "1000"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"gngs", "=", 
   RowBox[{"CompileGateNoise", "[", 
    RowBox[{"$gaussianQubitGateSet", ",", 
     RowBox[{"nmm", "[", 
      RowBox[{"0.9", ",", "100", ",", "1000"}], "]"}], ",", "1"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"protocol", "=", 
   RowBox[{"RBProtocol", "[", 
    RowBox[{"mlist", ",", "50", ",", "1", ",", 
     RowBox[{
      RowBox[{"TP", "[", "UI", "]"}], "/", "2"}], ",", 
     RowBox[{"0.99", 
      RowBox[{"TP", "[", "UI", "]"}]}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"survivals", "=", 
   RowBox[{"SimulateProtocol", "[", 
    RowBox[{"gngs", ",", "protocol", ",", 
     RowBox[{
     "SimulationExportName", "\[Rule]", "\"\<../data/nonmarkovian\>\""}]}], 
    "]"}]}], ";"}]}], "Input"],

Cell[BoxData["\<\"Saved ../data/nonmarkovian.h5\"\>"], "Print"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"model", "=", 
   RowBox[{
    RowBox[{
     RowBox[{"(", 
      RowBox[{"A", "-", "B"}], ")"}], 
     SuperscriptBox["p", "m"]}], "+", "B"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"normalizedSurvivals", "=", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"SequenceLengths", "[", "protocol", "]"}], ",", 
      RowBox[{
       RowBox[{"Total", "[", 
        RowBox[{
         RowBox[{"survivals", "[", 
          RowBox[{"[", 
           RowBox[{"All", ",", "1", ",", "All", ",", "1"}], "]"}], "]"}], ",", 
         RowBox[{"{", "2", "}"}]}], "]"}], "/", 
       RowBox[{
        RowBox[{"NumSequenceDraws", "[", "protocol", "]"}], "[", "]"}]}]}], 
     "}"}], "\[Transpose]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"fit", "=", 
   RowBox[{"FindFit", "[", 
    RowBox[{"normalizedSurvivals", ",", 
     RowBox[{"{", 
      RowBox[{"model", ",", 
       RowBox[{"0", "<", "p", "<", "1"}], ",", 
       RowBox[{"0", "<", "A", "<", "1"}], ",", 
       RowBox[{"0", "<", "B", "<", "1"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"A", ",", "B", ",", "p"}], "}"}], ",", "m"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"modelf", "=", 
   RowBox[{"Function", "[", 
    RowBox[{
     RowBox[{"{", "m", "}"}], ",", 
     RowBox[{"Evaluate", "[", 
      RowBox[{"model", "/.", "fit"}], "]"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Print", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"A", ",", "B", ",", "p"}], "}"}], "/.", "fit"}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Plot", "[", 
  RowBox[{
   RowBox[{"modelf", "[", "m", "]"}], ",", 
   RowBox[{"{", 
    RowBox[{"m", ",", "0", ",", 
     RowBox[{"Max", "@", 
      RowBox[{"SequenceLengths", "@", "protocol"}]}]}], "}"}], ",", 
   RowBox[{"Epilog", "\[Rule]", 
    RowBox[{"Map", "[", 
     RowBox[{"Point", ",", "normalizedSurvivals"}], "]"}]}], ",", 
   RowBox[{"PlotRange", "\[Rule]", 
    RowBox[{"{", 
     RowBox[{"All", ",", 
      RowBox[{"{", 
       RowBox[{"0", ",", "1"}], "}"}]}], "}"}]}]}], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0.9900172643598686`", ",", "0.5135724886618581`", ",", 
   "0.9957293874219489`"}], "}"}]], "Print"],

Cell[BoxData[
 GraphicsBox[{{{}, {}, 
    TagBox[
     {RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6], Opacity[
      1.], LineBox[CompressedData["
1:eJwVzWs0lAkcx/HZwWayFKVcdt3mJIklKhXm/1ciGqy246h1KUJZx7jVyExN
GmaEmSGVzbaZkfuW2lKhScXsKk6RdjvkcuaZQTcaBtVG2mdf/M73fF79bGNY
u+KoFAoliNz/faudZ/rFzzKKeNvPb2rQQFrVUXdf2TOwu/tAa16ngUyDiXhz
mRpM7kXofr6kAU5rNUVXNg06ptKgh+c18LmUH/SvVAezdBgH03M10KqrnZ6W
LkOvP/K8LSM1YOW5xn9KSkeTh1OsK4YaCH6umNRI3ZFpMPpEnv4OsJhvsD3D
F6NPPOdHeE6A8XTYteLYYBQwu6J3uI4D405m6LGlu7GfccyrFN6CNL31ymD4
HrQ7FNXo8cMbGA4Y6p3zi8IaTompkP8a4nZmTcR078eo5NvjEU2vwKflg253
8gGs5V52cKW+AoMmgmJJxKOiM+dEbPxLWDvln+yVdwgdKs8WSrrGYEuzxX1v
ahL2+Q9HO/uPwdWA6FV/sZJxmJf+qL97FO4uXrpgfCoFWcstPK33jIJFU/XA
tfZUHJANBJROj4Adfy61xSwdx2xvzSnEIxDlRjstTMhAbYC6PXz9CBivCapU
lB9Gnk2ans5LNbC8Q1hJ2iOY8cmEY12mBt/zQqec7zOxi73myL5wNShMjld8
nXsUf661rKevVAM7xbmC+U8W1s4zzKoHVSAzXXaGuoSLTb2cLmm9Cth9txj0
2GNo4dqzSJOpAoWqaOCG9DiG0YRPU4NV4PuUTtXv4+HejlCRvYMKmpJbUtRn
TmCJ3ICxhKqCU2bH+7PcstH0gUctlyDA9an0ov9QNtrv4t6xf0DA5Fy5w/uc
k/heFMNdXUlAUmDL7xG2fHRLHxzaW0DA3Rr3go+dfOwQ6zU2phHQ6fAwlnko
B7nCzc6ySAKmLDd/dVI3F3t2bGLLAgmY+I42N345F+8nzDCaPQio2jYzsSRQ
gNlM82/Eqwlwy2R4H5wRoFFi4GLeSgI4PZdyM0uEuMIrbqOCRsBQ2ejpRJc8
DKh50nHgsxJMMwTzxX152ODcaLhySgm01K3RqUdOocC9oaV4TAm93EgrgW0+
EoblCptBJdz8tLuC+ygfz9ntN1H0KuHHCwnWa+MLUBpDjWvrVEIHs/+j0dJC
nLrBbH7cpoT3NvFd9Q2F6DnKC6XIlQAJ3TEOgSLkhYrkwzeVECUvDRvQitBs
h9W6P68qYdMGzja/IjHqsZ70i+qV8Jyz1VGvRIz2Pp/y2aQPt9GMFWfFeFYr
9N1P+nrIL8M+v4px60/ynvWknRJvHvWuFmNw5LqgwTol2JRrrq6Xi7Gtq+iZ
E2ka7cC39FdiTLwvpPxdo4TakLU6qjdifDz+9nYrab9z2tfSCTE2pISw60jn
0LNvW02LcSBL34hHesH74i6LBTEW3nlX5kham9aXZ7xcgrHXJ1341Uooai5n
9ayQ4CrXkbhk0i6UhDCJuQTXsWsu7SGdJJqlG1pLML+yeKMr6bEak3v6jhK8
UFzRPFRF/r3rr+pwkmBwQpvzI9L0DbJCgYsEq+rK6hpJ72tz2au7QYKtLw42
FZBe0P+A7R4SZH8J38km/VtI6+qTWyTYPj87GkPa61yukY+3BHdPeuUFk34x
yJz9AhKc+ejotoX0f30DWkI=
       "]]},
     Annotation[#, "Charting`Private`Tag$2002995#1"]& ]}, {}, {}},
  AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
  Axes->{True, True},
  AxesLabel->{None, None},
  AxesOrigin->{0, 0},
  DisplayFunction->Identity,
  Epilog->{
    PointBox[{1, 0.9886522667788268}], 
    PointBox[{2, 0.9860971132468052}], 
    PointBox[{4, 0.9816171710139636}], 
    PointBox[{8, 0.9741926821136361}], 
    PointBox[{16, 0.9557418604252544}], 
    PointBox[{32, 0.9287617789122083}], 
    PointBox[{63, 0.8775952652859437}], 
    PointBox[{126, 0.8003332535057917}], 
    PointBox[{251, 0.6644342539174641}], 
    PointBox[{501, 0.5767919737590572}], 
    PointBox[{1000, 0.5178841693438324}]},
  Frame->{{False, False}, {False, False}},
  FrameLabel->{{None, None}, {None, None}},
  FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
  GridLines->{None, None},
  GridLinesStyle->Directive[
    GrayLevel[0.5, 0.4]],
  ImagePadding->All,
  Method->{
   "DefaultBoundaryStyle" -> Automatic, "DefaultMeshStyle" -> 
    AbsolutePointSize[6], "ScalingFunctions" -> None, 
    "CoordinatesToolOptions" -> {"DisplayFunction" -> ({
        (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
         Part[#, 1]], 
        (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
         Part[#, 2]]}& ), "CopiedValueFunction" -> ({
        (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
         Part[#, 1]], 
        (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
         Part[#, 2]]}& )}},
  PlotRange->{All, {0, 1}},
  PlotRangeClipping->True,
  PlotRangePadding->{{
     Scaled[0.02], 
     Scaled[0.02]}, {0, 0}},
  Ticks->{Automatic, Automatic}]], "Output"]
}, Open  ]]
}, Closed]]
}, Open  ]]
},
WindowSize->{1280, 1000},
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
CellTagsIndex->{
 "Info993710147481-3174582"->{
  Cell[153739, 3414, 2963, 42, 78, "Print",
   CellTags->"Info993710147481-3174582"]}
 }
*)
(*CellTagsIndex
CellTagsIndex->{
 {"Info993710147481-3174582", 231403, 5486}
 }
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1486, 35, 52, 0, 66, "Chapter"],
Cell[1541, 37, 184, 4, 33, "Text"],
Cell[1728, 43, 259, 7, 55, "Text"],
Cell[1990, 52, 279, 6, 78, "Input"],
Cell[2272, 60, 241, 7, 34, "Input"],
Cell[CellGroupData[{
Cell[2538, 71, 39, 0, 65, "Section"],
Cell[2580, 73, 158, 3, 30, "Text"],
Cell[2741, 78, 619, 19, 34, "Input"],
Cell[3363, 99, 1827, 49, 191, "Input"],
Cell[5193, 150, 95, 2, 33, "Text"],
Cell[5291, 154, 1253, 37, 101, "Input"],
Cell[CellGroupData[{
Cell[6569, 195, 271, 9, 34, "Input"],
Cell[6843, 206, 189, 6, 50, "Output"]
}, Open  ]]
}, Closed]],
Cell[CellGroupData[{
Cell[7081, 218, 43, 0, 51, "Section"],
Cell[CellGroupData[{
Cell[7149, 222, 2527, 68, 370, "Input"],
Cell[9679, 292, 148, 4, 34, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[9864, 301, 541, 16, 57, "Input"],
Cell[10408, 319, 79, 0, 24, "Print"]
}, Open  ]],
Cell[CellGroupData[{
Cell[10524, 324, 357, 10, 34, "Input"],
Cell[10884, 336, 9829, 251, 1375, "Output"]
}, Open  ]]
}, Closed]],
Cell[CellGroupData[{
Cell[20762, 593, 43, 0, 51, "Section"],
Cell[CellGroupData[{
Cell[20830, 597, 3151, 82, 462, "Input"],
Cell[23984, 681, 147, 4, 34, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24168, 690, 541, 16, 57, "Input"],
Cell[24712, 708, 79, 0, 24, "Print"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24828, 713, 357, 10, 34, "Input"],
Cell[25188, 725, 51165, 918, 1398, "Output"]
}, Open  ]]
}, Closed]],
Cell[CellGroupData[{
Cell[76402, 1649, 58, 0, 51, "Section"],
Cell[CellGroupData[{
Cell[76485, 1653, 3700, 100, 507, "Input"],
Cell[80188, 1755, 128, 3, 34, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[80353, 1763, 557, 16, 57, "Input"],
Cell[80913, 1781, 95, 1, 24, "Print"]
}, Open  ]],
Cell[CellGroupData[{
Cell[81045, 1787, 357, 10, 34, "Input"],
Cell[81405, 1799, 35679, 676, 1381, "Output"]
}, Open  ]]
}, Closed]],
Cell[CellGroupData[{
Cell[117133, 2481, 43, 0, 51, "Section"],
Cell[CellGroupData[{
Cell[117201, 2485, 4921, 144, 507, "Input"],
Cell[122125, 2631, 79, 0, 24, "Print"]
}, Open  ]],
Cell[CellGroupData[{
Cell[122241, 2636, 197, 6, 34, "Input"],
Cell[122441, 2644, 145, 4, 34, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[122623, 2653, 1101, 34, 79, "Input"],
Cell[123727, 2689, 13557, 282, 827, "Output"]
}, Open  ]]
}, Closed]],
Cell[CellGroupData[{
Cell[137333, 2977, 27, 0, 51, "Section"],
Cell[137363, 2979, 144, 3, 30, "Text"],
Cell[CellGroupData[{
Cell[137532, 2986, 1005, 27, 79, "Input"],
Cell[138540, 3015, 1087, 25, 57, "Output"],
Cell[139630, 3042, 841, 22, 48, "Output"]
}, Open  ]],
Cell[140486, 3067, 2353, 60, 192, "Input"],
Cell[CellGroupData[{
Cell[142864, 3131, 363, 10, 57, "Input"],
Cell[143230, 3143, 2853, 57, 119, "Output"]
}, Open  ]],
Cell[146098, 3203, 5092, 142, 442, "Input"],
Cell[CellGroupData[{
Cell[151215, 3349, 84, 1, 32, "Input"],
Cell[151302, 3352, 2343, 54, 103, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[153682, 3411, 54, 1, 32, "Input"],
Cell[153739, 3414, 2963, 42, 78, "Print",
 CellTags->"Info993710147481-3174582"]
}, Open  ]],
Cell[156717, 3459, 45, 0, 32, "Input"],
Cell[156765, 3461, 1317, 32, 239, "Input"],
Cell[CellGroupData[{
Cell[158107, 3497, 73, 1, 32, "Input"],
Cell[158183, 3500, 7054, 203, 286, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[165274, 3708, 1012, 25, 148, "Input"],
Cell[166289, 3735, 293, 7, 34, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[166619, 3747, 693, 18, 57, "Input"],
Cell[167315, 3767, 64, 0, 24, "Print"]
}, Open  ]],
Cell[CellGroupData[{
Cell[167416, 3772, 1964, 60, 149, "Input"],
Cell[169383, 3834, 142, 4, 26, "Print"],
Cell[169528, 3840, 3631, 77, 247, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[173196, 3922, 357, 10, 34, "Input"],
Cell[173556, 3934, 17924, 398, 1697, "Output"]
}, Open  ]]
}, Closed]],
Cell[CellGroupData[{
Cell[191529, 4338, 26, 0, 51, "Section"],
Cell[191558, 4340, 5672, 166, 552, "Input"],
Cell[197233, 4508, 110, 3, 33, "Text"],
Cell[197346, 4513, 2503, 75, 216, "Input"],
Cell[199852, 4590, 65, 0, 33, "Text"],
Cell[CellGroupData[{
Cell[199942, 4594, 3482, 101, 443, "Input"],
Cell[203427, 4697, 150, 4, 34, "Output"]
}, Open  ]],
Cell[203592, 4704, 4130, 103, 579, "Input"],
Cell[CellGroupData[{
Cell[207747, 4811, 1274, 35, 103, "Input"],
Cell[209024, 4848, 74, 0, 24, "Print"]
}, Open  ]],
Cell[CellGroupData[{
Cell[209135, 4853, 67, 1, 32, "Input"],
Cell[209205, 4856, 118, 3, 34, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[209360, 4864, 742, 21, 101, "Input"],
Cell[210105, 4887, 2216, 47, 239, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[212358, 4939, 377, 12, 34, "Input"],
Cell[212738, 4953, 672, 18, 66, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[213447, 4976, 251, 7, 55, "Input"],
Cell[213701, 4985, 31, 0, 32, "Output"],
Cell[213735, 4987, 31, 0, 32, "Output"]
}, Open  ]]
}, Closed]],
Cell[CellGroupData[{
Cell[213815, 4993, 31, 0, 51, "Section"],
Cell[213849, 4995, 809, 23, 73, "Input"],
Cell[CellGroupData[{
Cell[214683, 5022, 8687, 232, 832, "Input"],
Cell[223373, 5256, 276, 6, 34, "Output"]
}, Open  ]],
Cell[223664, 5265, 276, 9, 34, "Input"],
Cell[CellGroupData[{
Cell[223965, 5278, 31, 0, 32, "Input"],
Cell[223999, 5280, 182, 4, 34, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[224218, 5289, 799, 23, 80, "Input"],
Cell[225020, 5314, 63, 0, 24, "Print"]
}, Open  ]],
Cell[CellGroupData[{
Cell[225120, 5319, 2119, 64, 149, "Input"],
Cell[227242, 5385, 143, 4, 26, "Print"],
Cell[227388, 5391, 3576, 76, 250, "Output"]
}, Open  ]]
}, Closed]]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature mxpvoWOZcC@wjCwcYlFAeNdT *)
