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
    @Published var goodFuncToPlot : [plotDataType] = []
    @Published var energyFunctional : [plotDataType] = []
    @Published var potentialPlot : [plotDataType] = []
    
    /// eulerSolve:
    /// Solve the 1D Schrodinger Equation using euler's method
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
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
    
    /// rk4Solve:
    /// Solve the 1D Schrodinger Equation using rk4
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
    func rk4Solve(a: Double, steps: Double, Vf: PotentialList, ic: InitialCondition) {
        let precision = 1e-5
        
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
        
        for energyVal in stride(from: 0, to: 5, by: 0.2) {
//        let energyVal = 4.0 * 1.23
            var curPsiList = [ic.psi], curPsipList = [ic.psip]
            for i in 1..<xs.count {
                let V = vs[i], x = xs[i]
                
                let k0 = h * psiIter(psi: curPsi, psip: curPsip, xval: x),
                    l0 = h * psiPrimeIter(psi: curPsi, psip: curPsip, xval: x, energy: energyVal, V: V)
                
                let k1 = h * psiIter(psi: curPsi+(k0/2), psip: curPsip+(l0/2), xval: x+(h/2)),
                    l1 = h * psiPrimeIter(psi: curPsi+(k0/2), psip: curPsip+(l0/2), xval: x+(h/2), energy: energyVal, V: V)
                
                let k2 = h * psiIter(psi: curPsi+(k1/2), psip: curPsip+(l1/2), xval: x+(h/2)),
                    l2 = h * psiPrimeIter(psi: curPsi+(k1/2), psip: curPsip+(l1/2), xval: x+(h/2), energy: energyVal, V: V)
                
                let k3 = h * psiIter(psi: curPsi+h, psip: curPsip+h, xval: x+h),
                    l3 = h * psiPrimeIter(psi: curPsi+h, psip: curPsip+h, xval: x+h, energy: energyVal, V: V)
                
                let nextPsi = curPsi + (k0 + 2*k1 + 2*k2 + k3)/6,
                    nextPsiP = curPsip + (l0 + 2*l1 + 2*l2 + l3)/6
                
                curPsi  = nextPsi
                curPsip = nextPsiP
                
                curPsiList.append(curPsi)
                curPsipList.append(curPsip)
            }
            
            psiCollection.append(curPsiList)
            psipCollection.append(curPsipList)
            lastPointForBC.append((psi: curPsi, energy: energyVal))
            
            // check for sign change in psiLast

        }
        
//        print(curPsiList)
        
        // Boundary Condition is that psi(a) = 0 at the correct energy
        var energyFunc :[(psi: Double, energy: Double)] = []
        for i in 0..<lastPointForBC.count {
            let energyVal = lastPointForBC[i].energy
            let psiVal = lastPointForBC[i].psi
            if (abs(psiVal) < 0.1) {
                goodEnergyValues.append(energyVal)
                goodPsis.append(psiCollection[i])
            }
            energyFunc.append((psi: psiVal, energy: energyVal))
        }
        
        toPlotData(xvals: xs, yvals: psiCollection[6])
        fillPotentialPlot(potential: Vf)
        fillEnergyFunc(vals: energyFunc)
    }
    
//    func bisectionRootFinding(energyFunctional: [Double], leftGuess: Int, rightGuess: Int, tol: Double) -> Int{
//        var zeroVal = 1.0, nStep = 1
//        var leftVal = leftGuess, rightVal = rightGuess
//
//        assert(leftGuess > rightGuess)
//        while(abs(zeroVal) > tol) {
//            nStep += 1
//            let leftFunVal = energyFunctional[leftVal]
//            let rightFunVal = energyFunctional[rightVal]
//            let mid = (leftVal + rightVal) / 2
//            switch (sign(x : leftFunVal * zeroVal)) {
//            case 1:
//                leftVal = mid
//                break
//            case -1:
//                rightVal = mid
//                break
//            case 0:
//                rootidx = mid
//                break
//            default:
//                exit(123)
//            }
//        }
//    }
    
    /// sign:
    /// returns the sign of the function, 1, 0 or -1
    /// - Parameter x: Number to consider
    func sign(x: Double) -> Int {
        if(x == 0.0) {
            return 0
        } else if (x > 0.0) {
            return 1
        } else if (x < 0.0) {
            return -1
        }
        return 0
    }
    
    
    /// psiIter, psiPrimeIter
    /// Functions aiding the RK4 solving method, in order to make code more reuseable
    func psiIter(psi: Double, psip: Double, xval: Double) -> Double {
        return psip
    }
    
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
//        assert((xvals.count == yvals.count) && yvals.count > 0)
        for i in 0..<xvals.count {
            let x = xvals[i], y = yvals[i]
            goodFuncToPlot.append([.X: x, .Y: y])
        }
    }
    
    func toPlotData(xvals: [Double], yvals: [[Double]]) {
//        assert((xvals.count == yvals.count) && yvals.count > 0)
        print(yvals.count)
        for vals in yvals {
            for i in 0..<xvals.count {
                let x = xvals[i], y = vals[i]
                goodFuncToPlot.append([.X: x, .Y: y])
            }
        }
    }
    
    func fillEnergyFunc(vals: [(psi: Double, energy: Double)]) {
//        assert((xvals.count == yvals.count) && yvals.count > 0)
        for (psi, energy) in vals {
//            let x = tup.energy, y = tup.psi
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
