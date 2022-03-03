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
    
    var body: some View {
        HStack {
            Button("Solve", action: self.calculate)
                .frame(width: 100)
                .padding()

            CorePlot(dataForPlot: $solver.goodFuncToPlot,
                     changingPlotParameters: $plotData.plotArray[0].changingPlotParameters)
                .setPlotPadding(left: 10)
                .setPlotPadding(right: 10)
                .setPlotPadding(top: 10)
                .setPlotPadding(bottom: 10)
                .padding()
        }
    }
        
    
    func calculate() {
        let a = 2.0
        solver.rk4Solve(a: a, steps: 1000, Vf: squareWell(xMin: 0, xMax: a, steps: 1000, height: 0.0), ic: (psi: 0, psip: 1))
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
