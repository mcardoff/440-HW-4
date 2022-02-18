//
//  Potential.swift
//  Schrodinger-Solver
//
//  Created by Michael Cardiff on 2/18/22.
//

import Foundation

typealias CoordTuple = (x: Double, y: Double)
typealias ddFunc = (_: Double) -> Double

func generalWell(xMin: Double, xMax: Double, step: Double, f: ddFunc) -> [CoordTuple] {
    var ret : [CoordTuple] = []
    for xVal in stride(from: xMin+step, to: xMax-step, by: step) {
        ret.append((x: xVal, y: f(xVal)))
    }
    return ret
}

func squareWell(xMin: Double, xMax: Double, step: Double, height: Double) -> [CoordTuple] {
    return generalWell(xMin: xMin, xMax: xMax, step: step, f: {(_:Double) -> Double in return height})
}

func linearWell(xMin: Double, xMax: Double, step: Double, slope: Double) -> [CoordTuple] {
    return generalWell(xMin: xMin, xMax: xMax, step: step, f: {(x:Double) -> Double in return slope * x})
}

func quadraticWell(xMin: Double, xMax: Double, step: Double, amplitude: Double) -> [CoordTuple] {
    return generalWell(xMin: xMin, xMax: xMax, step: step, f: {(x:Double) -> Double in return amplitude * x * x})
}

func triangleWell(xMin: Double, xMax: Double, step: Double, amplitude: Double) -> [CoordTuple] {
    func triangleWell(x: Double) -> Double {
        if (x > 0.2*xMax && (x < xMax / 2.0)) {
            return amplitude * x
        } else if (x > xMax / 2.0 && x < 0.6*xMax) {
            return -amplitude * x + (amplitude * xMax / 2)
        } else {
            return 0.0
        }
    }
    return generalWell(xMin: xMin, xMax: xMax, step: step, f: triangleWell)
}
