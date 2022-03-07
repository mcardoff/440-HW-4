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
    @State var selector = 0
    
    // changed variables
    @State var numSteps: Int? = 250
    @State var wellWidth: Double? = 2.0
    @State var potentialString: String = ""
    @State var potentialVal: PotentialType = .square
    
    private var intFormatter: NumberFormatter = {
        let f = NumberFormatter()
        f.numberStyle = .decimal
        return f
    }()
    
    private var doubleFormatter: NumberFormatter = {
        let f = NumberFormatter()
        f.minimumSignificantDigits = 2
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
                    Text("Potential")
                    Picker("", selection: $potentialVal) {
                        ForEach(PotentialType.allCases) {
                            potential in Text(potential.toString())
                        }
                    }.frame(width: 150.0)
                }.padding()
                
                Button("Solve", action: self.calculate)
                    .frame(width: 100)
                    .padding()
                
                Button("iter", action: self.iteratesel)
                    .frame(width: 100)
                    .padding()
            }
            
            CorePlot(dataForPlot: $solver.totalFuncToPlot[selector],
                     changingPlotParameters: $plotData.plotArray[0].changingPlotParameters)
                .setPlotPadding(left: 10)
                .setPlotPadding(right: 10)
                .setPlotPadding(top: 10)
                .setPlotPadding(bottom: 10)
                .padding()
        }
    }
        
    
    func iteratesel() {
        if selector < $solver.totalFuncToPlot.count - 1 {
            selector += 1
        } else {
            selector = 0
        }
    }
    
    func calculate() {
        let ic : InitialCondition = (psi: 0, psip: 1)
        let V = getPotential(xMin: 0, xMax: wellWidth!, steps: numSteps!, choice: potentialVal)
        solver.boundaryValProblem(a: wellWidth!, steps: numSteps!, Vf: V, ic: ic)
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
