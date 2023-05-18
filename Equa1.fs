namespace Equation
open LVec

type MEqu (N:int, equ:MFunc list)=
    member val E = equ
    member val N = N
        // let mutable maxi = 0
        // for func in equ do maxi <- System.Math.Max(maxi, func.N)
        // maxi

    member this.Jacobi (pts:float list, ?dx) = 
        let dx = defaultArg dx 1e-12
        List.map(fun (func:MFunc) -> 
            List.map(fun i -> 
                MFunc.Partial(func.F, i, pts, dx)) 
                    [0..(func.N-1)]) this.E

    member this.Solve (point0:float list, 
        ?precision, ?maxLoop, ?showOpt:ShowSolveOpt) = 
        let precision, maxLoop, showOpt =
            defaultArg precision 1e-3,
            defaultArg maxLoop 300,
            defaultArg showOpt ShowSolveOpt.OnlyResult
        let mutable pts, counter = point0, 0
        while counter <= maxLoop do
            // f(x) = (f1, f2, f3)^T; x = (a, b, c)^T
            let fVals = List.map(fun (mfnc:MFunc) -> mfnc.F pts) this.E
            let matFx, matJacobi = [fVals] |> T, this.Jacobi(pts)
            if abs(Det(matJacobi)) < precision then counter <- maxLoop + 2
            else
                let deltaX = 
                    (MutiMat(Rev(matJacobi), matFx) |> T).[0]
                pts <- Subtract(pts, deltaX)
                match showOpt with
                | ShowSolveOpt.WithDetail -> 
                    printfn "L%i: F[]%A = %A" counter pts (this.E |> List.map(fun (mfnc:MFunc) -> mfnc.F pts))
                    printf "Jacobi =\t"
                    PrtMat matJacobi
                    printf "DeltaX =\t"
                    PrtMat [deltaX]
                | ShowSolveOpt.WithProcess -> printfn "L%i\t%A" counter pts
                | _ -> ()
                if Mudulus deltaX <= precision then counter <- maxLoop + 2
                counter <- counter + 1
        pts,
        List.map(fun (mfnc:MFunc) -> mfnc.F pts) this.E,
        if counter = maxLoop + 2 then true else false
