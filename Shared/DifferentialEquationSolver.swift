//
//  DifferentialEquationSolver.swift
//  Schrodinger-Solver
//
//  Created by Michael Cardiff on 2/18/22.
//

import Foundation

typealias Potential = (_:Double) -> Double

// solve Schrodinger Equation in 1D with potential
func eulerSolve(a: Double, stepSize: Double, V: Potential) {
//    let stepSize = 0.1
    
    var psi : [Double] = []
    var xs : [Double] = []
    for x in stride(from: 0, to: a, by: stepSize) {
        // do n steps
        
    }
}
