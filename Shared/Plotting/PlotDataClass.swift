//
//  PlotDataClass.swift
//  SwiftUICorePlotExample
//  Shared
//
//  Created by Jeff Terry on 12/16/20.
//

import Foundation
import SwiftUI
import CorePlot

class PlotDataClass: NSObject, ObservableObject {
    
    @MainActor @Published var plotData = [plotDataType]()
    @MainActor @Published var changingPlotParameters: ChangingPlotParameters = ChangingPlotParameters()
    @MainActor @Published var calculatedText = ""
    //In case you want to plot vs point number
    @MainActor @Published var pointNumber = 1.0
    
    @MainActor init(fromLine line: Bool) {
        //Must call super init before initializing plot
        super.init()
        
        //Intitialize the first plot
        self.plotBlank()
    }
    
    @MainActor func plotBlank() {
        zeroData()
        
        //set the Plot Parameters
        changingPlotParameters.yMax = 1.0
        changingPlotParameters.yMin = -1.0
        changingPlotParameters.xMax = 2.1
        changingPlotParameters.xMin = -0.1
        changingPlotParameters.xLabel = "x"
        changingPlotParameters.yLabel = "y"
        changingPlotParameters.lineColor = .blue()
        changingPlotParameters.title = "1-D Schrodinger"
        
    }
    
    @MainActor func zeroData() {
        plotData = []
        pointNumber = 1.0
    }
    
    @MainActor func appendData(dataPoint: [plotDataType]) {
        plotData.append(contentsOf: dataPoint)
        pointNumber += 1.0
    }
}
