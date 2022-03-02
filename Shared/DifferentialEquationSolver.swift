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
typealias IterationFun = (_:Double, _:Double, _:Double) -> Double

class SchrodingerSolver: NSObject, ObservableObject {
    
    @Published var psiVal : [Double] = []
    @Published var psipVal : [Double] = []
    @Published var xVal : [Double] = []
    @Published var goodFuncsToPlot : [[plotDataType]] = []
    @Published var energyFunctional : [plotDataType] = []
    
    // solve Schrodinger Equation in 1D with potential
    func eulerSolve(a: Double, steps: Double, Vf: PotentialList, ic: InitialCondition) {
        let precision = 1e-5
        
        //things to return
        let xs : [Double] = Vf.xs
        let vs : [Double] = Vf.Vs
        var psiCollection : [[Double]] = [], psipCollection : [[Double]] = []
        let stepSize = a / steps, h = stepSize
        
        // iteration values
//        let xstride = stride(from: stepSize, to: a, by: stepSize)
        var curPsi = ic.psi, curPsip = ic.psip
        
        // root finding
        var lastPointForBC : [(psi: Double, energy: Double)] = []
        var goodEnergyValues : [Double] = [], goodPsis : [[Double]] = []
        
        for energyVal in stride(from: 0, to: 10, by: 0.25) {
            var curPsiList = [ic.psi]
            var curPsipList = [ic.psip]
            for i in 1..<xs.count {
                // do n steps in total
                let V = vs[i], x = xs[i]
                let nextPsi  =  curPsi + h * curPsip
                let nextPsiP = curPsip + (h * Schrod(x: x, mass: 1.0, hbar: 1.0, energy: energyVal, V: V) * curPsi)
                curPsi  = nextPsi
                curPsip = nextPsiP
                
                curPsiList.append(curPsi)
                curPsipList.append(curPsip)
            }
            psiCollection.append(curPsiList)
            psipCollection.append(curPsipList)
            lastPointForBC.append((psi: curPsi, energy: energyVal))
        }
        
        for i in 0..<lastPointForBC.count {
            let energyVal = lastPointForBC[i].energy
            let psiVal = lastPointForBC[i].psi
            if (abs(psiVal) < precision) {
                goodEnergyValues.append(energyVal)
                goodPsis.append(psiCollection[i])
            }
        }
        
        psiVal.append(contentsOf: psiCollection[4])
        xVal.append(contentsOf: xs)
        toPlotData(xvals: xVal, yvals: psiCollection)
//        toDataCollection(xvals: xs, funVals: goodPsis)
//        toEnergyFunc(vals: lastPointForBC)
    }
    
    // solve Schrodinger Equation in 1D with potential
    func rk4Solve(a: Double, steps: Double, Vf: PotentialList, ic: InitialCondition) {
        let precision = 1e-5
        
        //things to return
        let xs : [Double] = Vf.xs
        let vs : [Double] = Vf.Vs
        var psiCollection : [[Double]] = [], psipCollection : [[Double]] = []
        let stepSize = a / steps, h = stepSize
        
        // iteration values
//        let xstride = stride(from: stepSize, to: a, by: stepSize)
        var curPsi = ic.psi, curPsip = ic.psip
        
        // root finding
        var lastPointForBC : [(psi: Double, energy: Double)] = []
        var goodEnergyValues : [Double] = [], goodPsis : [[Double]] = []
        
        for energyVal in stride(from: 0, to: 1, by: 0.25) {
            var curPsiList = [ic.psi], curPsipList = [ic.psip]
            for i in 0..<xs.count {
                let V = vs[i], x = xs[i]
//                let energyVal : Double =
                let k0 = h * psiIter(psi: curPsi, psip: curPsip, xval: x),
                    l0 = h * psiPrimeIter(psi: curPsi, psip: curPsip, xval: x, energy: energyVal, V: V),
                    k1 = h * psiIter(psi: curPsi+(k0/2), psip: curPsip+(l0/2), xval: x+(h/2)),
                    l1 = h * psiPrimeIter(psi: curPsi+(k0/2), psip: curPsip+(l0/2), xval: x+(h/2), energy: energyVal, V: V)
                
                let k2 = h * psiIter(psi: curPsi+(k1/2), psip: curPsip+(l1/2), xval: x+(h/2))
                let l2 = h * psiPrimeIter(psi: curPsi+(k1/2), psip: curPsip+(l1/2), xval: x+(h/2), energy: energyVal, V: V)
                
                let k3 = h * psiIter(psi: curPsi+h, psip: curPsip+h, xval: x+h)
                let l3 = h * psiPrimeIter(psi: curPsi+h, psip: curPsip+h, xval: x+h, energy: energyVal, V: V)
                
                let nextPsi = curPsi + (k0 + 2*k1 + 2*k2 + k3)/6
                let nextPsiP = curPsip + (l0 + 2*l1 + 2*l2 + l3)/6
                
                curPsi  = nextPsi
                curPsip = nextPsiP
                
                curPsiList.append(curPsi)
                curPsipList.append(curPsip)
            }
            
            psiCollection.append(curPsiList)
            psipCollection.append(curPsipList)
            lastPointForBC.append((psi: curPsi, energy: energyVal))
            
        }
        
        // Boundary Condition is that psi(a) = 0 at the correct energy
        for i in 0..<lastPointForBC.count {
            let energyVal = lastPointForBC[i].energy
            let psiVal = lastPointForBC[i].psi
            if (abs(psiVal) < precision) {
                goodEnergyValues.append(energyVal)
                goodPsis.append(psiCollection[i])
            }
        }
        
//        psiVal.append(contentsOf: psiCollection[0])
//        xVal.append(contentsOf: xs)
//        toPlotData(xvals: xVal, yvals: psiCollection)
//        toDataCollection(xvals: xs, funVals: goodPsis)
//        toEnergyFunc(vals: lastPointForBC)
    }
    
    func psiIter(psi: Double, psip: Double, xval: Double) -> Double {
        return psip
    }
    
    func psiPrimeIter(psi: Double, psip: Double, xval: Double, energy: Double, V: Double) -> Double {
        return Schrod(x: xval, mass: 1.0, hbar: 1.0, energy: energy, V: V) * psi
    }
    
    func Schrod(x: Double, mass: Double, hbar: Double, energy: Double, V: Double) -> Double {
        let consts = -2.0 * mass / (hbar * hbar)
        let rest = energy - V
        return consts * rest
    }
    
    
    // deprecated versions of the above functions for when V was a function not list
    func psiPrimeIterFunctionV(psi: Double, psip: Double, xval: Double, energy: Double, V: PotentialFunc) -> Double {
        return SchrodFunctionV(x: xval, mass: 1.0, hbar: 1.0, energy: energy, V: V) * psi
    }
    
    func SchrodFunctionV(x: Double, mass: Double, hbar: Double, energy: Double, V: PotentialFunc) -> Double {
        let consts = -2.0 * mass / (hbar * hbar)
        let rest = energy - V(x)
        return consts * rest
    }
    
    func toPlotData(xvals: [Double], yvals: [[Double]]) {
//        assert((xvals.count == yvals.count) && yvals.count > 0)
        print(yvals.count)
        for vals in yvals {
            var tempList : [plotDataType] = []
            for i in 0..<xvals.count {
                let x = xvals[i], y = vals[i]
                tempList.append([.X: x, .Y: y])
            }
            goodFuncsToPlot.append(tempList)
        }
    }
    
    func fillEnergyFunc(vals: [(psi: Double, energy: Double)]) {
//        assert((xvals.count == yvals.count) && yvals.count > 0)
        for tup in vals {
            let x = tup.energy, y = tup.psi
            energyFunctional.append([.X: x, .Y: y])
        }
    }
}
