namespace Equation

type FuncSolveOpt =
| FindMaximum = 1
| FindMinimum = 0

type ShowSolveOpt =
| OnlyResult = 0
| WithProcess = 1
| WithDetail = 2

type MMS = // Math Service
    static member Partial (fu:float list -> double, vari:int, pt:float list, ?delta:double) =
            let dv = defaultArg delta 1e-6
            let pt_ = List.mapi(fun i x -> 
                if i = vari then x + dv
                else x) pt
            // printfn $"Partial {vari} at {pt} = {(fu pt_ - fu pt) / dv}"
            (fu pt_ - fu pt) / dv