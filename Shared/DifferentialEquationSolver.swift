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

typealias PotentialFunc = (_:Double) -> Double
typealias PotentialList = (xs: [Double], Vs: [Double])
typealias InitialCondition = (psi: Double, psip: Double)
typealias PhaseSpacePt = (psi: Double, psip: Double)
typealias Iterfunctype = (_: Double, _: Double, _: Double, _: Double, _: Double, _: Double) -> PhaseSpacePt

class SchrodingerSolver: NSObject, ObservableObject {
    
    @Published var goodFuncToPlot : [plotDataType] = []
    @Published var totalFuncToPlot : [[plotDataType]] = [[[.X:0.0, .Y:0.0]]]
    @Published var energyFunctional : [plotDataType] = []
    @Published var potentialPlot : [plotDataType] = []
    
    /// rknSolve:
    /// Solve the 1D Schrodinger Equation using any function
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
    func rknSolve(a: Double, steps: Double, Vf: PotentialList, ic: InitialCondition, iterfunc: Iterfunctype) {
        let precision = 1e-3
        
        //things to return
        let xs : [Double] = Vf.xs
        let vs : [Double] = Vf.Vs
        let stepSize = a / steps, h = stepSize
        var psiCollection : [[Double]] = [], psipCollection : [[Double]] = []
        
        // iteration values
        var curPsi = ic.psi, curPsip = ic.psip
        
        // root finding
        var lastPointForBC : [(psi: Double, energy: Double)] = []
        var goodEnergyValues : [Double] = [], goodPsis : [[Double]] = []
        
        for energyVal in stride(from: 1.0, to: 1.5, by: 0.1) {
//        let energyVal = Double.pi*Double.pi / 8
            var curPsiList = [ic.psi], curPsipList = [ic.psip]
            for i in 1..<xs.count {
                let V = vs[i], x = xs[i]
                
                let deltatup : PhaseSpacePt = iterfunc(h, curPsi, curPsip, x, energyVal, V)
                
                let nextPsi = curPsi + deltatup.psi,
                    nextPsiP = curPsip + deltatup.psip
                
                curPsi  = nextPsi
                curPsip = nextPsiP
                
                curPsiList.append(curPsi)
                curPsipList.append(curPsip)
            }
            
            psiCollection.append(curPsiList)
            psipCollection.append(curPsipList)
            lastPointForBC.append((psi: curPsi, energy: energyVal))
            
            // reset the state dumbass
            curPsi = ic.psi
            curPsip = ic.psip
            
            // check for sign change in psiLast...
            // find two close vals that have a change in sign

        }

        var energyFunc : [(psi: Double, energy: Double)] = []
        for i in 0..<lastPointForBC.count { //
            let energyVal = lastPointForBC[i].energy
            let psiVal = lastPointForBC[i].psi
            if (abs(psiVal) < precision) {
                goodEnergyValues.append(energyVal)
                goodPsis.append(psiCollection[i])
            }
            
            
            
            energyFunc.append((psi: psiVal, energy: energyVal))
        }
        
        toPlotData(xvals: xs, yvals: psiCollection)
        fillPotentialPlot(potential: Vf)
        fillEnergyFunc(vals: energyFunc)
    }
    
    /// eulerSolve:
    /// Solve the 1D Schrodinger Equation using euler's method
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
    func eulerSolve(a: Double, steps: Double, Vf: PotentialList, ic: InitialCondition) {
        /// euler
        /// Matching iterfunc for the parent function
        func euler(h: Double, psi: Double, psip: Double, x: Double, energy: Double, V: Double) -> PhaseSpacePt {
            let k0 = h * psiIter(psi: psi, psip: psip, xval: x),
                l0 = h * psiPrimeIter(psi: psi, psip: psip, xval: x, energy: energy, V: V)
            
            return (psi: k0, psip: l0)
        }
        
        rknSolve(a: a, steps: steps, Vf: Vf, ic: ic, iterfunc: euler)
    }
    
    /// rk4Solve:
    /// Solve the 1D Schrodinger Equation using rk4
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
    func rk4Solve(a: Double, steps: Double, Vf: PotentialList, ic: InitialCondition) {
        /// rk4
        /// Matching iterfunc for the above function
        func rk4(h: Double, psi: Double, psip: Double, x: Double, energy: Double, V: Double) -> PhaseSpacePt {
            let k0 = h * psiIter(psi: psi, psip: psip, xval: x),
                l0 = h * psiPrimeIter(psi: psi, psip: psip, xval: x, energy: energy, V: V)
            
            let k1 = h * psiIter(psi: psi+(k0/2), psip: psip+(l0/2), xval: x+(h/2)),
                l1 = h * psiPrimeIter(psi: psi+(k0/2), psip: psip+(l0/2), xval: x+(h/2), energy: energy, V: V)
            
            let k2 = h * psiIter(psi: psi+(k1/2), psip: psip+(l1/2), xval: x+(h/2)),
                l2 = h * psiPrimeIter(psi: psi+(k1/2), psip: psip+(l1/2), xval: x+(h/2), energy: energy, V: V)
            
            let k3 = h * psiIter(psi: psi+h, psip: psip+h, xval: x+h),
                l3 = h * psiPrimeIter(psi: psi+h, psip: psip+h, xval: x+h, energy: energy, V: V)
            
            return (psi: (k0 + 2*k1 + 2*k2 + k3)/6, psip: (l0 + 2*l1 + 2*l2 + l3)/6)
        }
        
        rknSolve(a: a, steps: steps, Vf: Vf, ic: ic, iterfunc: rk4)
        
        if(self.totalFuncToPlot.isEmpty) {
            print("EMPTY")
        }
    }
    
    /// psiIter
    /// Functions aiding the RK4 solving method, in order to make code more reuseable
    func psiIter(psi: Double, psip: Double, xval: Double) -> Double {
        return psip
    }
    
    /// psiPrimeIter
    /// Functions aiding the RK4 solving method, in order to make code more reuseable
    func psiPrimeIter(psi: Double, psip: Double, xval: Double, energy: Double, V: Double) -> Double {
        return Schrod(x: xval, mass: 1.0, hbar: 1.0, energy: energy, V: V) * psi
    }
    
    /// Schrod
    /// Specific implementation of the schrodinger equation when turned into a system and solved for the derivative term:
    func Schrod(x: Double, mass: Double, hbar: Double, energy: Double, V: Double) -> Double {
        let consts = -2.0 * mass / (hbar * hbar)
        let rest = energy - V
        return consts * rest
    }
    
    /// Helper functions to turn data values in RK4/Euler to plottable things
    func toPlotData(xvals: [Double], yvals: [Double]) {
        for i in 0..<xvals.count {
            let x = xvals[i], y = yvals[i]
            goodFuncToPlot.append([.X: x, .Y: y])
        }
    }
    
    func toPlotData(xvals: [Double], yvals: [[Double]]) {
        totalFuncToPlot = []
        for vals in yvals {
            var tempList : [plotDataType] = []
            for (x,y) in zip(xvals, vals) {
                tempList.append([.X: x, .Y: y])
                goodFuncToPlot.append([.X: x, .Y: y])
            }
            totalFuncToPlot.append(tempList)
        }
    }
    
    func fillEnergyFunc(vals: [(psi: Double, energy: Double)]) {
        for (psi, energy) in vals {
            energyFunctional.append([.X: energy, .Y: psi])
        }
    }
    
    func fillPotentialPlot(potential: PotentialList) {
        for (x, V) in zip(potential.xs, potential.Vs) {
            potentialPlot.append([.X: x, .Y: V])
        }
    }
    
    
    // deprecated versions of the iteration functions for when V was a function not list
    func psiPrimeIterFunctionV(psi: Double, psip: Double, xval: Double, energy: Double, V: PotentialFunc) -> Double {
        return SchrodFunctionV(x: xval, mass: 1.0, hbar: 1.0, energy: energy, V: V) * psi
    }
    
    func SchrodFunctionV(x: Double, mass: Double, hbar: Double, energy: Double, V: PotentialFunc) -> Double {
        let consts = -2.0 * mass / (hbar * hbar)
        let rest = energy - V(x)
        return consts * rest
    }
}
