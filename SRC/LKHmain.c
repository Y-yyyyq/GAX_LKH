#include "LKH.h"
#include "Genetic.h"

/*
 * This file contains the main function of the program.
 */

int main(int argc, char *argv[])
{
    GainType Cost, OldOptimum; // Cost: 当前计算的旅行路线的成本。OldOptimum: 记录上一次的最优解。
    double Time, LastTime;     // 用于跟踪时间的变量

    /* Read the specification of the problem */
    if (argc >= 2) // 读取参数文件名: 如果命令行参数数量大于或等于2，则将第一个参数作为参数文件名。
        ParameterFileName = argv[1];
    ReadParameters(); // 读取参数: 调用函数从参数文件中读取配置，如问题类型、最大运行时间等。
    StartTime = LastTime = GetTime();
    MaxMatrixDimension = 20000; // 防止内存溢出
    MergeWithTour =
        Recombination == GPX2 ? MergeWithTourGPX2 : Recombination == CLARIST ? MergeWithTourCLARIST
                                                                             : MergeWithTourIPT;
    ReadProblem(); // 从文件或其他来源读取具体的旅行商问题数据，包括城市的位置和连接。

    if (SubproblemSize > 0)
    {
        if (DelaunayPartitioning)
            SolveDelaunaySubproblems();
        else if (KarpPartitioning)
            SolveKarpSubproblems();
        else if (KCenterPartitioning)
            SolveKCenterSubproblems();
        else if (KMeansPartitioning)
            SolveKMeansSubproblems();
        else if (RohePartitioning)
            SolveRoheSubproblems();
        else if (MoorePartitioning || SierpinskiPartitioning)
            SolveSFCSubproblems();
        else
            SolveTourSegmentSubproblems(); // 分割方法
        return EXIT_SUCCESS;
    }
    AllocateStructures(); // 为算法运行所需的数据结构分配内存
    CreateCandidateSet(); // 创建候选集: 准备一个候选集，通常用于优化算法，以便快速访问可能的解
    InitializeStatistics();

    if (Norm != 0)
        BestCost = PLUS_INFINITY; // 将最佳成本设为正无穷，表示尚未找到解
    else
    {
        /* The ascent has solved the problem! */
        Optimum = BestCost = (GainType)LowerBound; // 否则，将当前的最优解设置为下界
        UpdateStatistics(Optimum, GetTime() - LastTime);
        RecordBetterTour();
        RecordBestTour();
        WriteTour(OutputTourFileName, BestTour, BestCost);
        WriteTour(TourFileName, BestTour, BestCost);
        Runs = 0; // 重置运行次数: 如果找到解，重置运行次数。
    }

    /* Find a specified number (Runs) of local optima */
    for (Run = 1; Run <= Runs; Run++)
    {
        LastTime = GetTime();
        if (LastTime - StartTime >= TotalTimeLimit)
        {
            if (TraceLevel >= 1)
                printff("*** Time limit exceeded ***\n");
            Run--;
            break;
        }
        Cost = FindTour(); /* using the Lin-Kernighan heuristic */ // 使用Lin-Kernighan启发式算法找到当前的最优解
        if (MaxPopulationSize > 1)
        { // 如果种群大小大于1,使用遗传算法
            /* Genetic algorithm */
            int i;
            for (i = 0; i < PopulationSize; i++)
            { // 循环处理种群中的每个个体: 遍历整个种群。
                GainType OldCost = Cost;
                Cost = MergeTourWithIndividual(i);  // 合并当前路径与个体: 将当前路径与种群中的个体合并，更新成本。
                if (TraceLevel >= 1 && Cost < OldCost)
                { // 输出合并结果: 如果合并后的成本低于原来的成本，并根据调试级别输出合并信息。
                    printff("  Merged with %d: Cost = " GainFormat, i + 1,
                            Cost);
                    if (Optimum != MINUS_INFINITY && Optimum != 0) // 计算并打印当前成本与最优解之间的差距
                        printff(", Gap = %0.4f%%",
                                100.0 * (Cost - Optimum) / Optimum);
                    printff("\n");
                }
            }
            if (!HasFitness(Cost))
            { // 检查成本是否有适应度: 如果当前成本没有适应度，则进行进一步处理。
                if (PopulationSize < MaxPopulationSize)
                {
                    AddToPopulation(Cost); // 将解添加到种群: 如果种群未满，添加当前解到种群中。
                    if (TraceLevel >= 1)
                        PrintPopulation(); // 打印当前种群信息: 如果调试级别足够，打印种群状态。
                }
                else if (Cost < Fitness[PopulationSize - 1])
                {

                    i = ReplacementIndividual(Cost); // 如果种群已满，检查当前成本是否低于最差个体的适应度

                    ReplaceIndividualWithTour(i, Cost); // 找到适应度最差的个体并用当前路径替换。

                    if (TraceLevel >= 1)
                        PrintPopulation(); // 打印种群信息: 输出更新后的种群信息。
                }
            }
        }
        else if (Run > 1)
            Cost = MergeTourWithBestTour(); // 合并当前路径与最佳路径: 如果不是第一次运行，则将当前路径与历史最佳路径合并。
        if (Cost < BestCost)
        {
            BestCost = Cost; // 更新最佳成本: 如果当前成本低于历史最佳成本，进行更新。
            RecordBetterTour();
            RecordBestTour();
            WriteTour(OutputTourFileName, BestTour, BestCost);
            WriteTour(TourFileName, BestTour, BestCost);
        }
        OldOptimum = Optimum; // 记录上一次的最优解以便后续比较。
        if (Cost < Optimum)
        { // 如果当前成本低于已知的最优解。
            if (FirstNode->InputSuc)
            {
                Node *N = FirstNode;
                while ((N = N->InputSuc = N->Suc) != FirstNode)
                    ; // 如果存在输入后继节点，更新节点链表。
            }
            Optimum = Cost;
            printff("*** New optimum = " GainFormat " ***\n", Optimum);
        }
        Time = fabs(GetTime() - LastTime); // 计算当前运行的持续时间。
        UpdateStatistics(Cost, Time);
        if (TraceLevel >= 1 && Cost != PLUS_INFINITY)
        { // 输出当前运行的状态: 如果调试级别足够且成本不是无穷大，输出当前运行的信息
            printff("Run %d: Cost = " GainFormat, Run, Cost);
            if (Optimum != MINUS_INFINITY && Optimum != 0)
                printff(", Gap = %0.4f%%",
                        100.0 * (Cost - Optimum) / Optimum); // 输出与最优解的差距: 计算并打印当前成本与最优解之间的差距。
            printff(", Time = %0.2f sec. %s\n\n", Time,      // 输出当前运行的时间和状态（小于、等于或大于最优解）
                    Cost < Optimum ? "<" : Cost == Optimum ? "="
                                                           : "");
        }
        if (StopAtOptimum && Cost == OldOptimum && MaxPopulationSize >= 1)
        { // 如果设定了停止条件且当前成本等于旧的最优解，则准备结束运行。
            Runs = Run;
            break; // 调整运行次数并退出循环: 将运行次数设置为当前循环次数并退出
        }
        if (PopulationSize >= 2 &&
            (PopulationSize == MaxPopulationSize ||
             Run >= 2 * MaxPopulationSize) &&
            Run < Runs)
        { // 检查种群大小与运行次数: 如果种群大小大于或等于2，且达到种群的最大大小或已运行两倍的最大大小，则进行交叉操作。
            Node *N;
            int Parent1, Parent2;
            Parent1 = LinearSelection(PopulationSize, 1.25); // 选择父代: 使用线性选择方法选择父代个体。
            do
                Parent2 = LinearSelection(PopulationSize, 1.25);
            while (Parent2 == Parent1); // 确保不同父代: 确保选择的两个父代不同。
            ApplyCrossover(Parent1, Parent2); // 应用交叉操作: 将选择的两个父代进行交叉生成新个体
            N = FirstNode;                    // 遍历节点: 遍历所有节点，更新候选集
            do
            {
                if (ProblemType != HCP && ProblemType != HPP)
                { // 检查问题类型: 如果问题不是HCP（哈密尔顿回路）或HPP（哈密尔顿路径），则添加候选。
                    int d = C(N, N->Suc);
                    AddCandidate(N, N->Suc, d, INT_MAX);
                    AddCandidate(N->Suc, N, d, INT_MAX);
                } // 添加候选: 根据节点间的距离添加候选。
                N = N->InitialSuc = N->Suc;
            } while (N != FirstNode);
        } // 继续遍历节点直到回到起始节点。
        SRandom(++Seed); // 更新随机种子，以便在每次运行中引入随机性
    }
    PrintStatistics();
    system("pause");
    return EXIT_SUCCESS;
}
