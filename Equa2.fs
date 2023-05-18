namespace Equation
open LVec

type FuncSolveOpt =
| FindMaximum = 1
| FindMinimum = 0

type ShowSolveOpt =
| OnlyResult = 0
| WithProcess = 1
| WithDetail = 2

type MFunc(n:int, f:float list -> double) = class
    member val F = f
    member val N = n

    static member Partial (fu:float list -> double, vari:int, pt:float list, ?delta:double) =
        let dv = defaultArg delta 1e-8
        let pt_ = List.mapi(fun i x -> 
            if i = vari then x + dv
            else x) pt
        (fu pt_ - fu pt) / dv
    
    member this.Extremum (point0:float list, step0, ?precision0, ?maxLoops:int, ?findOpt:FuncSolveOpt, ?showOpt:ShowSolveOpt) =
        let findOpt, showOpt = defaultArg findOpt FuncSolveOpt.FindMinimum, defaultArg showOpt ShowSolveOpt.OnlyResult
        let Fun =
            match findOpt with
            | FuncSolveOpt.FindMaximum -> this.F >> (fun x -> -x)
            | _ -> this.F
        // Solve for Minimum below.
        let loops = defaultArg maxLoops 1000
        let precision = defaultArg precision0 1e-6
        let mutable pts = point0
        let mutable step = step0
        let mutable counter, flag = 0, true
        let mutable grad = List.map(fun id -> MFunc.Partial(Fun, id, pts, precision)) [0..(this.N - 1)]
        
        // Stop when step < d_x
        while flag && counter <= loops do
            let mutable inFlag = true
            while flag && inFlag do
                let dPt = Muti(-step, grad)
                if Mudulus dPt < precision then flag <- false
                let nxtPt, semiNxtPt = Add(pts, dPt), Add(pts, Muti(0.5, dPt))
                if Fun nxtPt < Fun pts && Fun nxtPt < Fun semiNxtPt then inFlag <- false
                else step <- step / 2.
            pts <- Add(pts, Muti(-step, grad))
            
            match showOpt with
            | ShowSolveOpt.WithProcess 
                -> printfn "L%i: %A %f" counter pts (this.F pts)
            | ShowSolveOpt.WithDetail
                -> printfn "L%i: F%A = %f, Grad%A, Step %f" counter pts (this.F pts) grad step
            | _ -> ()
            grad <- List.map(fun id -> MFunc.Partial(Fun, id, pts, precision)) [0..(this.N - 1)]
            counter <- counter + 1
        if flag then 
            // Failed - Break by the counter.
            pts, this.F pts, false
        else
            // Success - Break by the flag.
            pts, this.F pts, true

    member this.ExtremumWith(point0:float list, step0:double, ?equs:(float list -> double) list, ?ines:(float list -> double)list, 
        ?miu, ?rou, ?precision, ?maxLoops:int, ?findOpt:FuncSolveOpt, ?showOpt:ShowSolveOpt) =
        let rou, precision, maxLoops, findOpt, showOpt, equs, ines = 
            defaultArg rou 1.025,
            defaultArg precision 1e-4,
            defaultArg maxLoops 1000,
            defaultArg findOpt FuncSolveOpt.FindMinimum,
            defaultArg showOpt ShowSolveOpt.OnlyResult,
            defaultArg equs [fun (x:float list) -> 0.],
            defaultArg ines [fun (x:float list) -> infinity]
        let mutable miu1 = defaultArg miu 5.
        let mutable miu2 = miu1
        let func =
            match findOpt with
            | FuncSolveOpt.FindMaximum -> this.F >> (fun x -> -x)
            | _ -> this.F
        // Solve minimum below.
        let phi (x:float list):double =
            func x + 0.5 * miu1 * (List.map(fun f -> f x) equs |> (List.fold(fun acc elem -> acc + elem ** 2.)) 0.)
                 - miu2 * (List.map(fun f -> f x) ines |> (List.fold(fun acc elem -> acc + log elem)) 0.)
        
        let mutable counter, grad, pts, step = 0, [], point0, step0
        while counter <= maxLoops do
            grad <- List.map(fun id -> MFunc.Partial(phi, id, pts, precision)) [0..(this.N - 1)]
            let nP st = Add(pts, Muti(-st, grad))
            while not (phi(nP step) < phi (nP (step/2.)) && phi(nP step) < phi (pts)) 
                && Mudulus(Muti(step, grad)) >= precision do
                step <- step / 2.
            pts <- nP step
            match showOpt with
            | ShowSolveOpt.WithProcess
                -> printfn "L%i: %A %f" counter pts (this.F pts)
            | ShowSolveOpt.WithDetail 
                -> printfn "L%i: F|Φ%A = %.3f|%.3f, ▽%A, Step%.3f, μ %5f %f" counter pts (this.F pts) (phi pts) grad step miu1 miu2
            | _ -> ()
            miu1 <- miu1 * rou
            miu2 <- miu2 / rou
            step <- step * rou * rou
            counter <- counter + 1
            if Mudulus(Muti(step, grad)) < precision then counter <- maxLoops + 2
        pts, this.F pts, if counter = maxLoops + 2 then true else false

    override this.ToString() =
        $"MFunc(x[{this.N}])"
end
