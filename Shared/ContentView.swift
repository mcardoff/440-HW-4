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
    
    var body: some View {
        HStack {
            VStack {
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
        let a = 2.0
        let steps = 250
        let ic : InitialCondition = (psi: 0, psip: 1)
//        let V = squareWell(xMin: 0, xMax: a, steps: 1000, height: 0.0)
//        let V = linearWell(xMin: 0, xMax: a, steps: steps, slope: 14.0)
//        let V = quadraticWell(xMin: 0, xMax: a, steps: steps, amplitude: 1.0)
        let V = centeredQuadraticWell(xMin: 0, xMax: a, steps: steps, amplitude: 1.0)

        solver.boundaryValProblem(a: a, steps: steps, Vf: V, ic: ic)
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
