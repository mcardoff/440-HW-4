//
//  PlotClass.swift
//  Test Plot
//
//  Created by Jeff Terry on 1/25/21.
//

import Foundation

class PlotClass: ObservableObject {
    
    @Published var plotArray: [PlotDataClass]
    
    @MainActor init() {
        self.plotArray = [PlotDataClass.init(fromLine: true)]
        self.plotArray.append(contentsOf: [PlotDataClass.init(fromLine: true)])
        
    }
}
