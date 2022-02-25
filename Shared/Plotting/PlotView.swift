//
//  PlotView.swift
//  Shared
//
//  Created by Jeff Terry on 1/25/21.
//

import SwiftUI
import CorePlot

typealias plotDataType = [CPTScatterPlotField : Double]

struct PlotView: View {
    @EnvironmentObject var plotData :PlotClass
    
    @ObservedObject private var calculator = CalculatePlotData()
    @State var isChecked:Bool = false
    @State var tempInput = ""
    @State var selector = 0
    
    var body: some View {
        
        VStack{
            
            CorePlot(dataForPlot: $plotData.plotArray[selector].plotData, changingPlotParameters: $plotData.plotArray[selector].changingPlotParameters)
                .setPlotPadding(left: 10)
                .setPlotPadding(right: 10)
                .setPlotPadding(top: 10)
                .setPlotPadding(bottom: 10)
                .padding()
            
            Divider()
            
            HStack{
                Button("1s 1s", action: {
                    Task.init{
                        self.selector = 0
                        await self.calculateVsR(psil: "1s", psir: "1s")
                    }
                }).padding()
                
                Button("1s 2px", action: {
                    Task.init{
                        self.selector = 0
                        await self.calculateVsR(psil: "1s", psir: "2px")
                    }
                }).padding()
                
                Button("2px 2px", action: {
                    Task.init{
                        self.selector = 0
                        await self.calculateVsR(psil: "2px", psir: "2px")
                    }
                }).padding()
            }
        }
    }
    
    @MainActor func setupPlotDataModel(selector: Int){
        calculator.plotDataModel = self.plotData.plotArray[selector]
    }
    
    func calculateVsR(psil: String, psir: String) async {
        setupPlotDataModel(selector: 0)
        let _ = await withTaskGroup(of:  Void.self) { taskGroup in
            taskGroup.addTask {
//                var temp = 0.0
                //Calculate the new plotting data and place in the plotDataModel
                if psir == "1s" {
                    if psil == "2px" {
                        await calculator.plotOverlapVsR(psil: psi2px, psir: psi1s)
                    } else {
                        await calculator.plotOverlapVsR(psil: psi1s, psir: psi1s)
                    }
                } else if psir == "2px" {
                    if psil == "2px" {
                        await calculator.plotOverlapVsR(psil: psi2px, psir: psi2px)
                    } else {
                        await calculator.plotOverlapVsR(psil: psi1s, psir: psi2px)
                    }
                }
                
                // This forces a SwiftUI update. Force a SwiftUI update.
                await self.plotData.objectWillChange.send()
            }
        }
    }
    
    
    /// calculate
    /// Function accepts the command to start the calculation from the GUI
    func calculate() async {
        setupPlotDataModel(selector: 0)
        let _ = await withTaskGroup(of:  Void.self) { taskGroup in
            taskGroup.addTask {
//                var temp = 0.0
                //Calculate the new plotting data and place in the plotDataModel
                await calculator.ploteToTheMinusX()
                
                // This forces a SwiftUI update. Force a SwiftUI update.
                await self.plotData.objectWillChange.send()
            }
        }
    }
    
    /// calculate
    /// Function accepts the command to start the calculation from the GUI
    func calculate2() async {
        setupPlotDataModel(selector: 1)
        let _ = await withTaskGroup(of:  Void.self) { taskGroup in
            taskGroup.addTask {
//                var temp = 0.0
                //Calculate the new plotting data and place in the plotDataModel
                await calculator.plotYEqualsX()
                
                // This forces a SwiftUI update. Force a SwiftUI update.
                await self.plotData.objectWillChange.send()
            }
        }
    }
}
