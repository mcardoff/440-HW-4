//
//  DifferentialEquationSolver.swift
//  Schrodinger-Solver
//
//  Created by Michael Cardiff on 2/18/22.
//

/*
 d
 -- / x(t) \ = / f(x,y,t) \
 dt \ y(t) /   \ g(x,y,t) /
 
 for the Schrodinger equation:
 d(psi)
 ------ = psip
   dt

 d(psip)      m (E - V)
 ------- = -2 --------- psi
    dt          hb^2
 
 linearize derivative:
 
 psi(x+h) - psi(x)
 ----------------- = psip(x) => psi(x+h) = psi(x) + h (psip(x))
          h

 psip(x+h) - psip(x)      m (E - V)                                    /     m (E - V(x)) \
 ------------------- = -2 --------- psi(x) => psip(x+h) = psip(x) + h  | -2 ------------- |
           h                 hb^2                                      \       hb^2       /
 
 */

import Foundation

typealias Potential = (_:Double) -> Double
typealias InitialCondition = (psi: Double, psip: Double)

class SchrodingerSolver: NSObject, ObservableObject {
    
    @Published var psiVal : [Double] = []
    @Published var psipVal : [Double] = []
    @Published var xVal : [Double] = []
    @Published var zipped : [plotDataType] = []
    
    // solve Schrodinger Equation in 1D with potential
    func eulerSolve(a: Double, steps: Double, V: Potential, ic: InitialCondition) {
        //    let stepSize = 0.1
        
        var psi : [Double] = [ic.psi]
        var psip : [Double] = [ic.psip]
        var xs : [Double] = [0.0]
        let stepSize = a / steps
        var curPsi = ic.psi
        var curPsip = ic.psip
        for x in stride(from: stepSize, to: a-stepSize, by: stepSize) {
            // do n steps
            xs.append(x)
            let en = Double.pi * Double.pi / 8
            let nextPsi  =  curPsi + stepSize * curPsip
            let nextPsiP = curPsip + stepSize * Schrod(x: x, mass: 1.0, hbar: 1.0, energy: en, V: V) * curPsi
            curPsi  = nextPsi
            curPsip = nextPsiP
            
            psi.append(curPsi)
            psip.append(curPsip)
        }
        
        psi.append(0)
        xs.append(a)
        psiVal.append(contentsOf: psi)
        psipVal.append(contentsOf: psip)
        xVal.append(contentsOf: xs)
        toPlotData(xvals: xVal, yvals: psiVal)
    }
    
    func Schrod(x: Double, mass: Double, hbar: Double, energy: Double, V: Potential) -> Double {
        let consts = -2.0 * mass / (hbar * hbar)
        let rest = energy - V(x)
        return consts * rest
    }
    
    func toPlotData(xvals: [Double], yvals: [Double]) {
        assert((xvals.count == yvals.count) && yvals.count > 0)
        for i in 0..<xvals.count {
            let x = xvals[i], y = yvals[i]
            zipped.append([.X: x, .Y: y])
        }
    }
}
