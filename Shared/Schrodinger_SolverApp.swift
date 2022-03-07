//
//  Schrodinger_SolverApp.swift
//  Shared
//
//  Created by Michael Cardiff on 2/18/22.
//

import SwiftUI

@main
struct Schrodinger_SolverApp: App {
    
    @StateObject var plotData = PlotClass()
    
    var body: some Scene {
        WindowGroup {
            ContentView()
                .environmentObject(plotData)
        }
    }
}
