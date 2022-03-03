//
//  RootFinding.swift
//  Schrodinger-Solver
//
//  Created by Michael Cardiff on 3/2/22.
//

import Foundation

func bisectionRootFinding(energyFunctional: [Double], leftGuess: Int, rightGuess: Int, tol: Double) -> Int {
//    var zeroVal = 1.0, nStep = 1
//    var leftVal = leftGuess, rightVal = rightGuess
//
//    assert(leftGuess > rightGuess)
//    while(abs(zeroVal) > tol) {
//        nStep += 1
//        let leftFunVal = energyFunctional[leftVal]
//        let rightFunVal = energyFunctional[rightVal]
//        let mid = (leftVal + rightVal) / 2
//        switch (sign(x : leftFunVal * zeroVal)) {
//        case 1:
//            leftVal = mid
//            break
//        case -1:
//            rightVal = mid
//            break
//        case 0:
//            let rootidx = mid
//            break
//        default:
//            exit(123)
//        }
//    }
    return 0
}


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
