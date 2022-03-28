//
//  ContentView.swift
//  Shared
//
//  Created by Michael Cardiff on 2/18/22.
//

import SwiftUI
import CorePlot

typealias plotDataType = [CPTScatterPlotField : Double]

struct ContentView: View {
    
    @EnvironmentObject var plotData : PlotClass
    @ObservedObject var solver = SchrodingerSolver()
    @State var emptyPlot : [plotDataType] = []
    
    // changed variables
    @State var selector = 0
    @State var numSteps: Int? = 250
    @State var wellWidth: Double? = 2.0
    @State var amplitude: Double? = 0.5
    @State var eMin: Double? = 1.0
    @State var eMax: Double? = 30.0
    @State var eStep: Double? = 0.5
    @State var potentialVal: PotentialType = .square
    
    private var intFormatter: NumberFormatter = {
        let f = NumberFormatter()
        f.numberStyle = .decimal
        return f
    }()
    
    private var doubleFormatter: NumberFormatter = {
        let f = NumberFormatter()
        f.minimumSignificantDigits = 3
        f.maximumSignificantDigits = 9
        return f
    }()
    
    var body: some View {
        HStack {
            // variables to change: steps, a, potential, energy eigenval, energy search range
            VStack {
                VStack {
                    Text("Number of Steps")
                    TextField("Number of Steps in RK4", value: $numSteps, formatter: intFormatter)
                        .frame(width: 100.0)
                }.padding()
                
                VStack {
                    Text("Well Width")
                    TextField("Goes from [0,a]", value: $wellWidth, formatter: doubleFormatter)
                        .frame(width: 100.0)
                }.padding()
                
                VStack {
                    Text("Energy Min and Max")
                    HStack {
                        TextField("eMin", value: $eMin, formatter: doubleFormatter)
                            .frame(width: 100.0)
                        TextField("eMax", value: $eMax, formatter: doubleFormatter)
                            .frame(width: 100.0)
                    }
                }.padding()
                
                VStack {
                    Text("Energy Step")
                    HStack {
                        TextField("eStep", value: $eStep, formatter: doubleFormatter)
                            .frame(width: 100.0)
                    }
                }
                
                VStack {
                    Text("Potential")
                    Picker("", selection: $potentialVal) {
                        ForEach(PotentialType.allCases) {
                            potential in Text(potential.toString())
                        }
                    }.frame(width: 150.0)
                }.padding()
                
                VStack {
                    Text("Amplitude Parameter")
                    TextField("Amplitude Factor for Potential", value: $amplitude, formatter: doubleFormatter)
                        .frame(width: 100.0)
                }.padding()
                
                Button("Solve", action: self.calculate)
                    .frame(width: 100)
                    .padding()
                
                HStack {
                    Button("+E", action: self.increasesel)
                        .padding()
                    
                    Button("-E", action: self.decreasesel)
                        .padding()
                }
                Button("Clear", action: self.clear)
                    .frame(width: 100)
                    .padding()
            }
            TabView {
                CorePlot(
                    dataForPlot: $solver.totalFuncToPlot.count > 0 ? $solver.totalFuncToPlot[selector] : $emptyPlot,
                    changingPlotParameters: $plotData.plotArray[0].changingPlotParameters)
                    .setPlotPadding(left: 10)
                    .setPlotPadding(right: 10)
                    .setPlotPadding(top: 10)
                    .setPlotPadding(bottom: 10)
                    .padding()
                    .tabItem {
                        Text("Wavefunction Plot")
                    }
                
                CorePlot(
                    dataForPlot: $solver.potentialPlot,
                    changingPlotParameters: $plotData.plotArray[0].changingPlotParameters)
                    .setPlotPadding(left: 10)
                    .setPlotPadding(right: 10)
                    .setPlotPadding(top: 10)
                    .setPlotPadding(bottom: 10)
                    .padding()
                    .tabItem {
                        Text("Potential Plot")
                    }
                
                CorePlot(
                    dataForPlot: $solver.energyFunctional,
                    changingPlotParameters: $plotData.plotArray[0].changingPlotParameters)
                    .setPlotPadding(left: 10)
                    .setPlotPadding(right: 10)
                    .setPlotPadding(top: 10)
                    .setPlotPadding(bottom: 10)
                    .padding()
                    .tabItem {
                        Text("Energy Functional Plot")
                    }
                
                TextEditor(text: $solver.eigenvalueText)
                    .tabItem {
                        Text("Energy Eigenvalues")
                    }
            }
        }
    }
        
    func increasesel() {
        if selector < $solver.totalFuncToPlot.count - 1 {
            selector += 1
        } else {
            selector = 0
        }
    }
    
    func decreasesel() {
        if selector > 0 {
            selector -= 1
        } else {
            selector = 0
        }
    }
    
    func clear() {
        selector = 0
        solver.clearData()
    }
    
    func calculate() {
        self.clear()
        let ic : InitialCondition = (psi: 0, psip: 1)
        let V = getPotential(xMin: 0, xMax: wellWidth!, steps: numSteps!, choice: potentialVal, amplitude: amplitude!)
//        print(V)
        plotData.plotArray[0].changingPlotParameters.title = potentialVal.toString()
        solver.boundaryValProblem(a: wellWidth!, steps: numSteps!, Vf: V, ic: ic,
                                  eMin: eMin!, eMax: eMax!, eStride: eStep!)
        print(solver.totalFuncToPlot.count)
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
