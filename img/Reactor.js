function Reactor(optc,oppc,delta)
{
    var dr = 2,
    vrtotal = Math.PI*(dr**2)*dr/2,
    areaj = 2*Math.PI*dr**2,
    vr = 0.9*vrtotal,
    ca = 8.01,
    tj = 350,
    ch2 = 0.0,
    areahx = areaj,
    vj = 1*areaj,
    ccha = 0.0,
    vrca = vr*ca,
    vrch2 = vr*ch2,
    tr = 350,
    vrtr = vr*tr,
    vrccha = vr*ccha,
    pr;
    function ReactorMain( optc,oppc,delta)
    {
        var outPut = new Object();
        var k0 = 40000,  // preexponential factor 指数前系数(见公式2.10)
        e = 46.49e6, // activation energy
        t0 = 300,    // 输入氢气温度300K
        u = 851,     // 传热效率
        lambda = -190e6,  // 反应热
        roe = 801,      // 反应堆的密度
        ma = 93,        // 苯胺相对分子质量
        mcha = 99,      // 环己胺相对分子质量
        mh2 = 2,        // 氢气的相对分子质量
        cp = 3137,      // 比热容
        cj = 4183,      // 夹套内浓度
        roej = 1000,    // 夹套内的密度
        dr = 2,         // 反应器底径
        areatotal = areaj,             // 夹套内总面积
        tcin = 300,      // 冷却液的温度(K)
            // h2和cw的最大流量限
        fjmax = 0.005,// flow in cu m / sec
        fh2max = 0.035;// flow is kmol / sec
        /***************主要部分*********************** */
            // 模型主要部分
        // Integration  loop
        // 若反应器体积超过最大值，函数返回，直到反应器体积恢复过来
        if (vr > vrtotal )                                 // vr
        {
            outPut.TR = tr;
            outPut.PR = pr;
            return outPut;
        }    
            // 若苯胺明显要消耗殆尽，反应结束
        if (ca < 8.01 * 0.001)                               // ca
        {
            outPut.TR = tr;
            outPut.PR = pr;
            return outPut;
        }    
            // 由温度和成分计算压力
        var xa,xh2,xcha,psh2,psa,pscha,fj,fh2,q,dvr,dvrca,dvrch2,dtj,dvrtr,dvrccha
        xa = ca / (ca + ch2 + ccha);            // 反应器内苯胺占比
        xh2 = ch2 / (ca + ch2 + ccha);          // 反应器内氢气占比
        xcha = ccha / (ca + ch2 + ccha);        // 反应器内环己胺占比
        psh2 = 7300;                      // 氢气组分压力
        psa = Math.exp(11.6606 - 5329 / tr);       // 苯胺组分压力
        pscha = Math.exp(11.9125 - 4959 / tr);     // 环己胺组分压力
        pr = xa * psa + xh2 * psh2 + xcha * pscha;  // 反应器内总压力
        fj = fjmax * (1 - optc) * 2;            // 冷却液流量，由optc控制，其中(1 - optc) * 2相当于阀门开度
        fh2 = fh2max * oppc;                // 氢气输入流量，由oppc控制
        if (vr > vrtotal)
        {
            fh2 = 0;                      // 若反应器体积超过最大值，不再输入氢气
        }
        areahx = vr * areatotal / vrtotal;    // 当前冷却液有效接触面积
        q = u * areahx * (tr - tj);             // 当前冷却液吸收热量
        k = k0 * Math.exp(-e / tr / 8314);           // 产率系数(见公式2.10)
            // Derivative evaluations
        dvr = fh2 * mh2 / roe;                                       // 公式2.11(总质量平衡m3 / s)
        dvrca = -vr * k * ca * ch2;                                    // 公式2.12(苯胺的组分平衡kmol / s)
        dvrch2 = fh2 - 3 * vr * k * ca * ch2;                               // 公式2.13(氢的组分平衡kmol / s)
        dtj = fj * tcin / vj - fj * tj / vj + q / (cj * roej * vj);                 // 公式2.15(冷却液能量平衡J / s)
        dvrtr = fh2 * t0 * mh2 / roe - lambda * vr * k * ca * ch2 / roe / cp - q / (cp * roe);// 公式2.14(反应堆能量平衡J / s)
            // 更新当前反应器温度值与压力值，同时存储几个当前状态参数值，为下一周期做准备
        // Integration
        dvrccha = vr * k * ca * ch2;
        vr = vr + dvr * delta;               // 更新反应器体积
        vrca = vrca + dvrca * delta;
        vrch2 = vrch2 + dvrch2 * delta;
        vrccha = vrccha + dvrccha * delta;
        vrtr = vrtr + dvrtr * delta;
        tj = tj + dtj * delta;              // 更新夹套内温度
        tr = vrtr / vr;                 // 更新反应器温度
        ca = vrca / vr;                   // 更新反应器内苯胺浓度
        ch2 = vrch2 / vr;                 // 更新反应器内氢浓度
        ccha = vrccha / vr;               // 更新反应器内环己胺浓度
            // 输出当前反应器内温度及压力
        outPut.TR = tr;
        outPut.PR = pr;
        return outPut;
    }
    return ReactorMain;
}