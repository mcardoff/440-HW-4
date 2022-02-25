//
//  ChangingPlotParameters.swift
//  SwiftUICorePlotExample
//
//  Created by Jeff Terry on 12/19/20.
//

import SwiftUI
import CorePlot

class ChangingPlotParameters: NSObject, ObservableObject {
    
    //These plot parameters are adjustable
    
    var xLabel: String = "x"
    var yLabel: String = "y"
    var xMax : Double = 20.0
    var yMax : Double = 10.0
    var yMin : Double = 0.0
    var xMin : Double = 0.0
    var lineColor: CPTColor = .blue()
    var title: String = "Plot Title"
    
}
