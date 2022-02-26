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
    var xMax : Double = 2.1
    var yMax : Double = 1.0
    var yMin : Double = -1.0
    var xMin : Double = -0.1
    var lineColor: CPTColor = .blue()
    var title: String = "Plot Title"
    
}
