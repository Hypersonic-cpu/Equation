namespace Equation
open LVec

type MEqu (N:int, equ:(float list -> double) list)=
    member val E = equ
    member val N = N
        // let mutable maxi = 0
        // for func in equ do maxi <- System.Math.Max(maxi, func.N)
        // maxi

    member this.Jacobi (pts:float list, ?dx) = 
        let dx = defaultArg dx 1e-6 // Small dx will cause serious issues ~ 10^5
        List.map(fun (func:float list -> double) -> 
            List.map(fun i -> 
                MMS.Partial(func, i, pts, dx)) 
                    [0..(this.N-1)]) this.E

    member this.Solve (point0:float list, 
        ?precision, ?maxLoop, ?showOpt:ShowSolveOpt) = 
        let precision, maxLoop, showOpt =
            defaultArg precision 1e-4,
            defaultArg maxLoop 300,
            defaultArg showOpt ShowSolveOpt.OnlyResult
        let mutable pts, counter = point0, 0
        while counter <= maxLoop do
            // f(x) = (f1, f2, f3)^T; x = (a, b, c)^T
            let fVals = List.map(fun (mfnc:float list -> double) -> mfnc pts) this.E
            let matFx = [fVals] |> T
            let matJacobi = this.Jacobi(pts)
            PrtMat [pts]
            PrtMat matFx
            printf "Jacobi =\t"
            PrtMat matJacobi
            if abs(Det(matJacobi)) < precision then counter <- maxLoop + 1
            else
                let deltaX = 
                    (MutiMat(Rev(matJacobi), matFx) |> T).[0]
                pts <- Subtract(pts, deltaX)
                match showOpt with
                | ShowSolveOpt.WithDetail -> 
                    printfn "L%i: F[]%A = %A" counter pts (this.E |> List.map(fun (mfnc:float list -> double) -> mfnc pts))
                    printf "Jacobi =\t"
                    PrtMat matJacobi
                    printf "DeltaX =\t"
                    PrtMat [deltaX]
                | ShowSolveOpt.WithProcess -> printfn "L%i\t%A" counter pts
                | _ -> ()
                counter <- counter + 1
                if Mudulus deltaX <= precision then counter <- maxLoop + 2
        pts,
        List.map(fun (mfnc:float list -> double) -> mfnc pts) this.E,
        if counter = maxLoop + 2 then true else false
