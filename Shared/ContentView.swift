//
//  ContentView.swift
//  Shared
//
//  Created by Michael Cardiff on 2/18/22.
//

import SwiftUI

struct ContentView: View {
    
    @ObservedObject var solver = SchrodingerSolver()
    
    var body: some View {
        Button("Do Stuff", action: self.calculate)
            .frame(width: 100)
            .padding()
    }
    
    func calculate() {
        solver.eulerSolve(a: 10, steps: 100, V: {(_:Double) -> Double in return 0}, ic: (psi: 0, psip: 1))
        for item in solver.psiVal {
            print(item)
        }
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
