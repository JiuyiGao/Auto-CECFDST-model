# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
import math
from caeModules import *
from driverUtils import executeOnCaeStartup
import os

def parametric_model(fyso, fysi, fyl, fyv, fcc, fco, Do, Di, to, ti, dl, dv, H, B):
    # 直接使用传入的14个参数，转换为浮点数
    fyso = float(fyso)  # 外钢管屈服强度 (MPa)
    fysi = float(fysi)  # 内钢管屈服强度 (MPa)
    fyl = float(fyl)    # 纵筋屈服强度 (MPa)
    fyv = float(fyv)    # 箍筋屈服强度 (MPa)
    fcc = float(fcc)    # 夹层混凝土立方体抗压强度 (MPa)
    fco = float(fco)    # 外部混凝土立方体抗压强度 (MPa)
    Do = float(Do)      # 外钢管外径 (mm)
    Di = float(Di)      # 内钢管外径 (mm)
    to = float(to)      # 外钢管壁厚 (mm)
    ti = float(ti)      # 内钢管壁厚 (mm)
    dl = float(dl)      # 纵筋直径 (mm)
    dv = float(dv)      # 箍筋直径 (mm)
    H = float(H)        # 试件高度 (mm)
    B = float(B)        # 试件宽度 (mm)

    # 衍生参数
    Es1 = 206000.0                  # 外钢管弹性模量 (MPa)
    Es2 = 206000.0                  # 内钢管弹性模量 (MPa)
    Esl = 206000.0                  # 纵筋弹性模量 (MPa)
    Esv = 206000.0                  # 箍筋弹性模量 (MPa)
    us1 = 0.3                       # 外钢管泊松比
    us2 = 0.3                       # 内钢管泊松比
    usl = 0.3                       # 纵筋泊松比
    usv = 0.3                       # 箍筋泊松比
    fc1 = 0.79 * fcc                # 夹层混凝土棱柱体抗压强度 (MPa)
    fck1 = 0.6471 * fcc             # 夹层混凝土特征抗压强度 (MPa)
    fc2 = 0.79 * fco                # 外部混凝土棱柱体抗压强度 (MPa)
    Ec1 = 4730 * fc1 ** 0.5         # 夹层混凝土弹性模量 (MPa)
    Ec2 = 1e5 / (2.2 + 34.7 / fco)  # 外部混凝土弹性模量 (MPa)
    uc1 = 0.2                       # 夹层混凝土泊松比
    uc2 = 0.2                       # 外部混凝土泊松比
    Bo = B + B / 5                  # 端板宽度 (mm)
    d = H / 20                      # 端板厚度 (mm)
    U = H / 75                      # 施加位移 (mm)
    w = H / 30                      # 竖向网格

    # 外钢管本构
    Ee1 = 0.8 * fyso / Es1
    Enom_s1 = [1.5 * Ee1, 15 * Ee1, 150 * Ee1, 450 * Ee1]
    sigmanom_s1 = [fyso, fyso, 1.6 * fyso, 1.6 * fyso]
    Eture_s1 = [0] * 4
    sigmatrue_s1 = [0] * 4
    Ep_s1 = [0] * 4
    for i in range(4):
        Eture_s1[i] = math.log(1 + Enom_s1[i])
        sigmatrue_s1[i] = sigmanom_s1[i] * (1 + Enom_s1[i])
        Ep_s1[i] = Eture_s1[i] - sigmatrue_s1[i] / Es1
    Ep_s1[0] = 0
    COM_s1 = list(zip(sigmatrue_s1, Ep_s1))

    # 内钢管本构
    Ee2 = 0.8 * fysi / Es2
    Enom_s2 = [1.5 * Ee2, 15 * Ee2, 150 * Ee2, 450 * Ee2]
    sigmanom_s2 = [fysi, fysi, 1.6 * fysi, 1.6 * fysi]
    Eture_s2 = [0] * 4
    sigmatrue_s2 = [0] * 4
    Ep_s2 = [0] * 4
    for i in range(4):
        Eture_s2[i] = math.log(1 + Enom_s2[i])
        sigmatrue_s2[i] = sigmanom_s2[i] * (1 + Enom_s2[i])
        Ep_s2[i] = Eture_s2[i] - sigmatrue_s2[i] / Es2
    Ep_s2[0] = 0
    COM_s2 = list(zip(sigmatrue_s2, Ep_s2))

    # 夹层混凝土本构
    As = math.pi / 4 * (Do ** 2 - (Do - 2 * to) ** 2)
    Ac = math.pi / 4 * (Do - 2 * to) ** 2
    kexi1 = (As * fyso) / (Ac * fck1)
    e_c = (1300 + 12.5 * fc1) * 1e-6
    e_0 = e_c + 800 * kexi1 ** 0.2 * 1e-6
    p_0 = (2.36e-5) ** (0.25 + (kexi1 - 0.5) ** 7) * fc1 ** 0.5 * 0.5
    p_0 = max(p_0, 0.12)
    x1 = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
    length1 = len(x1)
    Enom_c1 = [0.0] * length1
    sigmanom_c1 = [0.0] * length1
    Etrue_c1 = [0.0] * length1
    sigmatrue_c1 = [0.0] * length1
    Ep_c1 = [0.0] * length1
    try:
        for i in range(length1):
            Enom_c1[i] = x1[i] * e_0
            x1_i = x1[i]
            if x1_i <= 1.0:
                sigmanom_c1[i] = fc1 * (2 * x1_i - x1_i ** 2)
            else:
                sigmanom_c1[i] = fc1 * x1_i / (p_0 * (x1_i - 1) ** 2 + x1_i)
            Etrue_c1[i] = math.log(1 + Enom_c1[i])
            sigmatrue_c1[i] = sigmanom_c1[i] * (1 + Enom_c1[i])
            Ep_c1[i] = Etrue_c1[i] - sigmatrue_c1[i] / Ec1
    except Exception as e:
        print('约束混凝土受压应力-应变计算错误:', str(e))
        raise
    start_idx = None
    for i in range(length1):
        if Ep_c1[i] > 0:
            start_idx = i
            break
    Enom_c1 = Enom_c1[start_idx:]
    sigmanom_c1 = sigmanom_c1[start_idx:]
    Etrue_c1 = Etrue_c1[start_idx:]
    sigmatrue_c1 = sigmatrue_c1[start_idx:]
    Ep_c1 = Ep_c1[start_idx:]
    Ep_c1[0] = 0
    COM_c1 = list(zip(sigmatrue_c1, Ep_c1))

    # 夹层混凝土断裂能
    f_cm = fck1 + 8
    d_max = 16.0
    f_t = 0.395 * fc1 ** 0.55
    G_F0 = 0.030
    f_cm0 = 10.0
    G_F = G_F0 * (f_cm / f_cm0) ** 0.7

    # 纵筋本构
    E_t1 = 0.01 * Esl
    e_y1 = fyl / Esl
    e_max1 = 0.1
    Enom_s3 = [e_y1, e_max1]
    sigmanom_s3 = [fyl, fyl + E_t1 * (e_max1 - e_y1)]
    length2 = len(Enom_s3)
    Etrue_s3 = [0.0] * length2
    sigmatrue_s3 = [0.0] * length2
    Ep_s3 = [0.0] * length2
    try:
        for i in range(length2):
            Etrue_s3[i] = math.log(1 + Enom_s3[i])
            sigmatrue_s3[i] = sigmanom_s3[i] * (1 + Enom_s3[i])
            Ep_s3[i] = Etrue_s3[i] - sigmatrue_s3[i] / Esl
        Ep_s3[0] = 0
    except Exception as e:
        print('Error in steel stress-strain calculation:', str(e))
        raise
    COM_s3 = list(zip(sigmatrue_s3, Ep_s3))

    # 箍筋本构
    E_t2 = 0.01 * Esv  # 强化阶段斜率 (206000 MPa)
    e_y2 = fyv / Esv  # 屈服应变
    e_max2 = 0.1  # 假设最大名义应变
    Enom_s4 = [e_y2, e_max2]  # 名义应变：屈服、强化
    sigmanom_s4 = [fyv, fyv + E_t2 * (e_max2 - e_y2)]  # 名义应力：屈服、强化
    length3 = len(Enom_s4)
    Etrue_s4 = [0.0] * length3
    sigmatrue_s4 = [0.0] * length3
    Ep_s4 = [0.0] * length3
    try:
        for i in range(length3):
            Etrue_s4[i] = math.log(1 + Enom_s4[i])  # 真实应变
            sigmatrue_s4[i] = sigmanom_s4[i] * (1 + Enom_s4[i])  # 真实应力
            Ep_s4[i] = Etrue_s4[i] - sigmatrue_s4[i] / Esv  # 塑性应变
        Ep_s4[0] = 0  # 屈服点塑性应变为 0
    except Exception as e:
        print('Error in stirrup stress-strain calculation:', str(e))
        raise
    COM_s4 = list(zip(sigmatrue_s4, Ep_s4))

    # 普通混凝土受压本构
    e_02 = (700 + 172 * math.sqrt(fc2)) * 1e-6  # 峰值应变
    n = Ec2 * e_02 / (Ec2 * e_02 - fc2)  # 上升段参数
    a_c = 0.157 * fc2 ** 0.785 - 0.905  # 下降段参数
    x2 = [0.2, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
    length4 = len(x2)
    Enom_c2 = [0.0] * length4
    sigmanom_c2 = [0.0] * length4
    Etrue_c2 = [0.0] * length4
    sigmatrue_c2 = [0.0] * length4
    Ep_c2 = [0.0] * length4
    try:
        for i in range(length4):
            Enom_c2[i] = x2[i] * e_02  # 名义应变
            x2_i = x2[i]
            if x2_i <= 1.0:
                # 上升阶段
                sigmanom_c2[i] = fc2 * (n * x2_i) / (n - 1 + x2_i ** n)
            else:
                # 下降阶段
                sigmanom_c2[i] = fc2 * x2_i / (a_c * (x2_i - 1) ** 2 + x2_i)
            Etrue_c2[i] = math.log(1 + Enom_c2[i])  # 真实应变
            sigmatrue_c2[i] = sigmanom_c2[i] * (1 + Enom_c2[i])  # 真实应力
            Ep_c2[i] = Etrue_c2[i] - sigmatrue_c2[i] / Ec2  # 塑性应变
    except Exception as e:
        print('Error in concrete compression stress-strain calculation:', str(e))
        raise
    start_idx = None
    for i in range(length4):
        if Ep_c2[i] > 0:
            start_idx = i
            break
    Enom_c2 = Enom_c2[start_idx:]
    sigmanom_c2 = sigmanom_c2[start_idx:]
    Etrue_c2 = Etrue_c2[start_idx:]
    sigmatrue_c2 = sigmatrue_c2[start_idx:]
    Ep_c2 = Ep_c2[start_idx:]
    Ep_c2[0] = 0  # 将第一个点的塑性应变设为 0
    COM_c2 = list(zip(sigmatrue_c2, Ep_c2))

    # 普通混凝土受拉本构
    f_tr = 0.395 * fc2 ** 0.55
    e_tr = 65 * f_tr ** 0.54 * 1e-6
    a_t = 0.312 * f_tr ** 2
    x3 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
    length5 = len(x3)
    Enom_t = [0.0] * length5
    sigmanom_t = [0.0] * length5
    Etrue_t = [0.0] * length5
    sigmatrue_t = [0.0] * length5
    Ep_t = [0.0] * length5
    try:
        for i in range(length5):
            Enom_t[i] = x3[i] * e_tr
            x3_i = x3[i]
            if x3_i <= 1.0:
                sigmanom_t[i] = f_tr
            else:
                sigmanom_t[i] = f_tr * x3_i / (a_t * (x3_i - 1) ** 1.7 + x3_i)
            Etrue_t[i] = math.log(1 + Enom_t[i])
            sigmatrue_t[i] = sigmanom_t[i] * (1 + Enom_t[i])
            Ep_t[i] = Etrue_t[i] - sigmatrue_t[i] / Ec2
        Ep_t[0] = 0
    except Exception as e:
        print('Error in concrete tension stress-strain calculation:', str(e))
        raise
    COM_t = list(zip(sigmatrue_t, Ep_t))

    # 内钢管建模
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(Di / 2, 0.0))
    p = mdb.models['Model-1'].Part(name='ngg', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseShellExtrude(sketch=s, depth=H)

    # 夹层混凝土建模
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(Di / 2, 0.0))
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(Do / 2, 0.0))
    p = mdb.models['Model-1'].Part(name='jchnt', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s1, depth=H)

    # 外钢管建模
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(Do / 2, 0.0))
    p = mdb.models['Model-1'].Part(name='wgg', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseShellExtrude(sketch=s, depth=H)

    # 外混凝土建模
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(Do / 2, 0.0))
    s1.rectangle(point1=(B / 2, B / 2), point2=(-B / 2, -B / 2))
    p = mdb.models['Model-1'].Part(name='whnt', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s1, depth=H)

    # 端板建模
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
    s.rectangle(point1=(Bo / 2, Bo / 2), point2=(-Bo / 2, -Bo / 2))
    p = mdb.models['Model-1'].Part(name='db', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=d)

    # 箍筋建模
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
    s1.rectangle(point1=((B - 50) / 2, (B - 50) / 2), point2=(-(B - 50) / 2, -(B - 50) / 2))
    p = mdb.models['Model-1'].Part(name='gj', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseWire(sketch=s1)

    # 纵筋建模
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
    s.Line(point1=(0.0, 0.0), point2=(H, 0.0))
    p = mdb.models['Model-1'].Part(name='zj', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseWire(sketch=s)

    # 网格划分
    # 纵筋网格
    p = mdb.models['Model-1'].parts['zj']
    p.seedPart(size=w, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    # 外混凝土网格
    p = mdb.models['Model-1'].parts['whnt']
    p.seedPart(size=w, deviationFactor=0.1, minSizeFactor=0.1)
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    v1, e, d1 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
        edge=e[12], rule=CENTER), point2=p.InterestingPoint(edge=e[10],
        rule=MIDDLE), point3=p.InterestingPoint(edge=e[11], rule=MIDDLE))
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#3 ]', ), )
    v2, e1, d2 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
        edge=e1[11], rule=CENTER), point2=p.InterestingPoint(edge=e1[23],
        rule=MIDDLE), point3=p.InterestingPoint(edge=e1[19], rule=MIDDLE))
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#f ]', ), )
    v1, e, d1 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
        edge=e[16], rule=MIDDLE), point2=p.InterestingPoint(edge=e[27],
        rule=MIDDLE), point3=p.InterestingPoint(edge=e[41], rule=MIDDLE))
    p.DatumPointByCoordinate(coords=((B - 50) / 2, (B - 50) / 2, 0.0))
    p.DatumPointByCoordinate(coords=((B - 50) / 2, -(B - 50) / 2, 0.0))
    p.DatumPointByCoordinate(coords=(-(B - 50) / 2, (B - 50) / 2, 0.0))
    p.DatumPointByCoordinate(coords=(-(B - 50) / 2, -(B - 50) / 2, 0.0))
    p.DatumPointByCoordinate(coords=((B - 50) / 2, (B - 50) / 2, H))
    p.DatumPointByCoordinate(coords=((B - 50) / 2, -(B - 50) / 2, H))
    p.DatumPointByCoordinate(coords=(-(B - 50) / 2, (B - 50) / 2, H))
    p.DatumPointByCoordinate(coords=(-(B - 50) / 2, -(B - 50) / 2, H))
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#d2 ]', ), )
    v2, e1, d2 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(point1=d2[13], point2=d2[12], point3=d2[8],
        cells=pickedCells)
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#e80 ]', ), )
    v1, e, d1 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(point1=d1[10], point2=d1[11], point3=d1[7],
        cells=pickedCells)
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#1b36 ]', ), )
    v2, e1, d2 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(point1=d2[12], point2=d2[10], point3=d2[8],
        cells=pickedCells)
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#fe8000 ]', ), )
    v1, e, d1 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(point1=d1[13], point2=d1[11], point3=d1[7],
        cells=pickedCells)
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=(
        '[#4a082402 #24a8000 #41480048 #10200050 #1162030 #2000410 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=7, constraint=FINER)
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=(
        '[#0 #1000008 #0 #40001000 #0 #9800200 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=8, constraint=FINER)
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=(
        '[#0 #800014 #0 #20002800 #80 #100 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=4, constraint=FINER)
    p.generateMesh()

    # 外钢管网格
    p = mdb.models['Model-1'].parts['wgg']
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=H / 2)
    p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    d2 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d2[2], faces=pickedFaces)
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#3 ]', ), )
    d1 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d1[3], faces=pickedFaces)
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#f ]', ), )
    d2 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d2[4], faces=pickedFaces)
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#f8888 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=8, constraint=FINER)
    p.seedPart(size=w, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    # 内钢管网格
    p = mdb.models['Model-1'].parts['ngg']
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=H / 2)
    p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    d1 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d1[2], faces=pickedFaces)
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#3 ]', ), )
    d2 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d2[3], faces=pickedFaces)
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#f ]', ), )
    d1 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d1[4], faces=pickedFaces)
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#f8888 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=8, constraint=FINER)
    p.seedPart(size=w, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    # 夹层混凝土网格
    p = mdb.models['Model-1'].parts['jchnt']
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    v2, e1, d2 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(point2=v2[0], point3=v2[1],
        cells=pickedCells, point1=p.InterestingPoint(edge=e1[0], rule=CENTER))
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#3 ]', ), )
    v1, e, d1 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
        edge=e[10], rule=CENTER), point2=p.InterestingPoint(edge=e[10],
        rule=MIDDLE), point3=p.InterestingPoint(edge=e[9], rule=MIDDLE))
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#f ]', ), )
    v2, e1, d2 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
        edge=e1[0], rule=MIDDLE), point2=p.InterestingPoint(edge=e1[26],
        rule=MIDDLE), point3=p.InterestingPoint(edge=e1[15], rule=MIDDLE))
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#45212000 #ce714 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=8, constraint=FINER)
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#10080000 #18e1 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=4, constraint=FINER)
    p.seedPart(size=w, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    # 箍筋网格
    p = mdb.models['Model-1'].parts['gj']
    p.seedPart(size=w, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    # 端板网格
    p = mdb.models['Model-1'].parts['db']
    p.seedPart(size=w, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    # 定义材料属性
    mdb.models['Model-1'].Material(name='db')
    mdb.models['Model-1'].materials['db'].Elastic(table=((100000000000000.0, 1e-14), ))

    mdb.models['Model-1'].Material(name='gj')
    mdb.models['Model-1'].materials['gj'].Elastic(table=((Esv, usv), ))
    mdb.models['Model-1'].materials['gj'].Plastic(scaleStress=None, table=(COM_s4))

    mdb.models['Model-1'].Material(name='zj')
    mdb.models['Model-1'].materials['zj'].Elastic(table=((Esl, usl), ))
    mdb.models['Model-1'].materials['zj'].Plastic(scaleStress=None, table=(COM_s3))

    mdb.models['Model-1'].Material(name='wgg')
    mdb.models['Model-1'].materials['wgg'].Elastic(table=((Es1, us1), ))
    mdb.models['Model-1'].materials['wgg'].Plastic(scaleStress=None, table=(COM_s1))

    mdb.models['Model-1'].Material(name='ngg')
    mdb.models['Model-1'].materials['ngg'].Elastic(table=((Es2, us2), ))
    mdb.models['Model-1'].materials['ngg'].Plastic(scaleStress=None, table=(COM_s2))

    mdb.models['Model-1'].Material(name='whnt')
    mdb.models['Model-1'].materials['whnt'].Elastic(table=((Ec2, uc2), ))
    mdb.models['Model-1'].materials['whnt'].ConcreteDamagedPlasticity(table=((30.0, 0.1, 1.16, 0.667, 0.0005), ))
    mdb.models['Model-1'].materials['whnt'].concreteDamagedPlasticity.ConcreteCompressionHardening(table=(COM_c2))
    mdb.models['Model-1'].materials['whnt'].concreteDamagedPlasticity.ConcreteTensionStiffening(table=(COM_t))

    mdb.models['Model-1'].Material(name='jchnt')
    mdb.models['Model-1'].materials['jchnt'].Elastic(table=((Ec1, uc1), ))
    mdb.models['Model-1'].materials['jchnt'].ConcreteDamagedPlasticity(table=((30.0, 0.1, 1.16, 0.667, 0.0005), ))
    mdb.models['Model-1'].materials['jchnt'].concreteDamagedPlasticity.ConcreteCompressionHardening(table=(COM_c1))
    mdb.models['Model-1'].materials['jchnt'].concreteDamagedPlasticity.ConcreteTensionStiffening(table=((f_t, G_F), ), type=GFI)

    # 创建截面
    mdb.models['Model-1'].HomogeneousSolidSection(name='db', material='db', thickness=None)
    mdb.models['Model-1'].HomogeneousShellSection(name='wgg', preIntegrate=OFF, material='wgg', thicknessType=UNIFORM, thickness=to, thicknessField='', nodalThicknessField='', idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, integrationRule=SIMPSON, numIntPts=9)
    mdb.models['Model-1'].HomogeneousShellSection(name='ngg', preIntegrate=OFF, material='ngg', thicknessType=UNIFORM, thickness=ti, thicknessField='', nodalThicknessField='', idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, integrationRule=SIMPSON, numIntPts=9)
    mdb.models['Model-1'].HomogeneousSolidSection(name='whnt', material='whnt', thickness=None)
    mdb.models['Model-1'].HomogeneousSolidSection(name='jchnt', material='jchnt', thickness=None)
    Al = (3.14 * dl ** 2) / 4
    mdb.models['Model-1'].TrussSection(name='zj', material='zj', area=Al)
    Av = (3.14 * dv ** 2) / 4
    mdb.models['Model-1'].TrussSection(name='gj', material='gj', area=Av)

    # 指派截面
    p = mdb.models['Model-1'].parts['db']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    region = regionToolset.Region(cells=cells)
    p.SectionAssignment(region=region, sectionName='db', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    p = mdb.models['Model-1'].parts['gj']
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#f ]', ), )
    region = regionToolset.Region(edges=edges)
    p.SectionAssignment(region=region, sectionName='gj', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    p = mdb.models['Model-1'].parts['jchnt']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#ff ]', ), )
    region = regionToolset.Region(cells=cells)
    p.SectionAssignment(region=region, sectionName='jchnt', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    p = mdb.models['Model-1'].parts['ngg']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#ff ]', ), )
    region = regionToolset.Region(faces=faces)
    p.SectionAssignment(region=region, sectionName='ngg', offset=0.0, offsetType=TOP_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    p = mdb.models['Model-1'].parts['wgg']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#ff ]', ), )
    region = regionToolset.Region(faces=faces)
    p.SectionAssignment(region=region, sectionName='wgg', offset=0.0, offsetType=TOP_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    p = mdb.models['Model-1'].parts['whnt']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#ffffffff ]', ), )
    region = regionToolset.Region(cells=cells)
    p.SectionAssignment(region=region, sectionName='whnt', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    p = mdb.models['Model-1'].parts['zj']
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#1 ]', ), )
    region = regionToolset.Region(edges=edges)
    p.SectionAssignment(region=region, sectionName='zj', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    # 装配
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-1'].parts['db']
    a.Instance(name='db-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['gj']
    a.Instance(name='gj-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['jchnt']
    a.Instance(name='jchnt-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['ngg']
    a.Instance(name='ngg-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['wgg']
    a.Instance(name='wgg-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['whnt']
    a.Instance(name='whnt-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['zj']
    a.Instance(name='zj-1', part=p, dependent=ON)
    a.translate(instanceList=('db-1', ), vector=(0.0, 0.0, H))
    p = mdb.models['Model-1'].parts['db']
    a.Instance(name='db-2', part=p, dependent=ON)
    a.translate(instanceList=('db-2', ), vector=(0.0, 0.0, -d))
    a.translate(instanceList=('gj-1', ), vector=(0.0, 0.0, 150.0))
    nu = int(math.ceil(H / 150)) - 1
    a.LinearInstancePattern(instanceList=('gj-1', ), direction1=(0.0, 0.0, 1.0),
        direction2=(0.0, 1.0, 0.0), number1=nu, number2=1, spacing1=150.0, spacing2=950.0)
    a.rotate(instanceList=('zj-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 10.0, 0.0), angle=-90.0)
    a.translate(instanceList=('zj-1', ), vector=(-(B - 50) / 2, -(B - 50) / 2, 0.0))
    a.LinearInstancePattern(instanceList=('zj-1', ), direction1=(1.0, 0.0, 0.0),
        direction2=(0.0, 1.0, 0.0), number1=5, number2=5, spacing1=(B - 50) / 4, spacing2=(B - 50) / 4)
    a.deleteFeatures(('zj-1-lin-2-2', 'zj-1-lin-2-3', 'zj-1-lin-2-4', ))
    a.deleteFeatures(('zj-1-lin-3-2', 'zj-1-lin-3-3', 'zj-1-lin-3-4', ))
    a.deleteFeatures(('zj-1-lin-4-2', 'zj-1-lin-4-3', 'zj-1-lin-4-4', ))

    # 钢筋笼合并
    a1 = mdb.models['Model-1'].rootAssembly
    all_instances = a1.instances.keys()
    target_instances = [inst for inst in all_instances if inst.startswith('gj') or inst.startswith('zj')]
    if not target_instances:
        raise ValueError("No instances starting with 'gj' or 'zj' found in the assembly.")
    a1.InstanceFromBooleanMerge(name='gjl', instances=tuple(a1.instances[inst] for inst in target_instances),
        originalInstances=SUPPRESS, domain=GEOMETRY)
    p = mdb.models['Model-1'].parts['gjl']
    p.seedPart(size=w, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    # 分析步
    a.regenerate()
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', maxNumInc=10000, initialInc=0.01, maxInc=0.1)

    # 创建相互作用属性
    mdb.models['Model-1'].ContactProperty('IntProp-1')
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((0.6, ), ),
        shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005, elasticSlipStiffness=None)
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)

    # 面面接触
    a1 = mdb.models['Model-1'].rootAssembly
    s1 = a1.instances['wgg-1'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#ff ]', ), )
    region1 = regionToolset.Region(side1Faces=side1Faces1)
    s1 = a1.instances['whnt-1'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#0:3 #884e00 #21 ]', ), )
    region2 = regionToolset.Region(side1Faces=side1Faces1)
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='wgg-whnt', createStepName='Step-1',
        main=region1, secondary=region2, sliding=FINITE, thickness=OFF, interactionProperty='IntProp-1',
        adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

    s1 = a1.instances['wgg-1'].faces
    side2Faces1 = s1.getSequenceFromMask(mask=('[#ff ]', ), )
    region1 = regionToolset.Region(side2Faces=side2Faces1)
    s1 = a1.instances['jchnt-1'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#21102830 #1 ]', ), )
    region2 = regionToolset.Region(side1Faces=side1Faces1)
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='wgg-jchnt', createStepName='Step-1',
        main=region1, secondary=region2, sliding=FINITE, thickness=OFF, interactionProperty='IntProp-1',
        adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

    a1.features['jchnt-1'].suppress()
    a1.features['jchnt-1'].resume()
    s1 = a1.instances['ngg-1'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#ff ]', ), )
    region1 = regionToolset.Region(side1Faces=side1Faces1)
    s1 = a1.instances['jchnt-1'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#80884380 #2 ]', ), )
    region2 = regionToolset.Region(side1Faces=side1Faces1)
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='ngg-jchnt', createStepName='Step-1',
        main=region1, secondary=region2, sliding=FINITE, thickness=OFF, interactionProperty='IntProp-1',
        adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

    # 创建约束
    a1 = mdb.models['Model-1'].rootAssembly
    e1 = a1.instances['gjl-1'].edges
    a1.Set(name='gjl-edges', edges=e1)
    edges1 = a1.sets['gjl-edges'].edges
    region1 = regionToolset.Region(edges=edges1)
    try:
        c1 = a1.instances['whnt-1'].cells
        if 'whnt-cells' in a1.sets:
            del a1.sets['whnt-cells']
        a1.Set(name='whnt-cells', cells=c1)
        host_region = regionToolset.Region(cells=a1.sets['whnt-cells'].cells)
    except KeyError:
        raise KeyError("Instance 'whnt-1' not found or has no cells. Please check the instance name or geometry.")
    mdb.models['Model-1'].EmbeddedRegion(name='gjl', embeddedRegion=region1, hostRegion=host_region,
        weightFactorTolerance=1e-06, absoluteTolerance=0.0, fractionalTolerance=0.05, toleranceMethod=BOTH)

    # 创建端板表面
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['db-1'].faces
    side1Faces1 = s1.findAt(((0, 0, H), ))
    a.Surface(side1Faces=side1Faces1, name='db-n-1')
    s1 = a.instances['db-2'].faces
    side1Faces1 = s1.findAt(((0, 0, 0), ))
    a.Surface(side1Faces=side1Faces1, name='db-n-2')
    s1 = a.instances['db-1'].faces
    side1Faces1 = s1.findAt(((0, 0, H + d), ))
    a.Surface(side1Faces=side1Faces1, name='db-w-1')
    s1 = a.instances['db-2'].faces
    side1Faces1 = s1.findAt(((0, 0, -d), ))
    a.Surface(side1Faces=side1Faces1, name='db-w-2')

    # 绑定约束
    region1 = a.surfaces['db-n-1']
    s1 = a.instances['whnt-1'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#a020800 #140480 #120090 #20600000 #40 ]', ), )
    s2 = a.instances['jchnt-1'].faces
    side1Faces2 = s2.getSequenceFromMask(mask=('[#42040000 #4 ]', ), )
    region2 = a.Surface(side1Faces=side1Faces1 + side1Faces2, name='hnt-1')
    mdb.models['Model-1'].Tie(name='db-hnt-1', main=region1, secondary=region2,
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    region1 = a.surfaces['db-n-2']
    s1 = a.instances['whnt-1'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#109200 #4601000 #840003 #82040000 #80 ]', ), )
    s2 = a.instances['jchnt-1'].faces
    side1Faces2 = s2.getSequenceFromMask(mask=('[#10600000 #8 ]', ), )
    region2 = a.Surface(side1Faces=side1Faces1 + side1Faces2, name='hnt-2')
    mdb.models['Model-1'].Tie(name='db-hnt-2', main=region1, secondary=region2,
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    # 壳-实体耦合约束
    s1 = a.instances['wgg-1'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#90880 ]', ), )
    s2 = a.instances['ngg-1'].edges
    side1Edges2 = s2.getSequenceFromMask(mask=('[#90880 ]', ), )
    region1 = a.Surface(side1Edges=side1Edges1 + side1Edges2, name='gg-1')
    region2 = a.surfaces['db-n-1']
    mdb.models['Model-1'].ShellSolidCoupling(name='db-gg-1', shellEdge=region1, solidFace=region2,
        positionToleranceMethod=COMPUTED)

    s1 = a.instances['wgg-1'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#68008 ]', ), )
    s2 = a.instances['ngg-1'].edges
    side1Edges2 = s2.getSequenceFromMask(mask=('[#68008 ]', ), )
    region1 = a.Surface(side1Edges=side1Edges1 + side1Edges2, name='gg-2')
    region2 = a.surfaces['db-n-2']
    mdb.models['Model-1'].ShellSolidCoupling(name='db-gg-2', shellEdge=region1, solidFace=region2,
        positionToleranceMethod=COMPUTED)

    # 点-面耦合
    a.ReferencePoint(point=(0.0, 0.0, -2 * d))
    a.ReferencePoint(point=(0.0, 0.0, H + 2 * d))
    r1 = a.referencePoints
    refPoints1 = (r1.findAt((0.0, 0.0, H + 2 * d)), )
    region1 = a.Set(referencePoints=refPoints1, name='m_Set-rp1')
    region2 = a.surfaces['db-w-1']
    mdb.models['Model-1'].Coupling(name='rp-1', controlPoint=region1, surface=region2,
        influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None,
        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

    r1 = a.referencePoints
    refPoints1 = (r1.findAt((0.0, 0.0, -2 * d)), )
    region1 = a.Set(referencePoints=refPoints1, name='m_Set-rp2')
    region2 = a.surfaces['db-w-2']
    mdb.models['Model-1'].Coupling(name='rp-2', controlPoint=region1, surface=region2,
        influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None,
        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

    # 位移施加
    r1 = a.referencePoints
    refPoints1 = (r1.findAt((0.0, 0.0, H + 2 * d)), )
    region = a.Set(referencePoints=refPoints1, name='Set-rp1')
    mdb.models['Model-1'].DisplacementBC(name='rp-1', createStepName='Initial',
        region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET,
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    r1 = a.referencePoints
    refPoints1 = (r1.findAt((0.0, 0.0, -2 * d)), )
    region = a.Set(referencePoints=refPoints1, name='Set-rp2')
    mdb.models['Model-1'].DisplacementBC(name='rp-2', createStepName='Initial',
        region=region, u1=SET, u2=SET, u3=UNSET, ur1=SET, ur2=SET, ur3=SET,
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    mdb.models['Model-1'].boundaryConditions['rp-2'].setValuesInStep(stepName='Step-1', u3=U)

    # 网格-钢筋笼指派单元类型-T3D2
    p = mdb.models['Model-1'].parts['gjl']
    elemType1 = mesh.ElemType(elemCode=T3D2, elemLibrary=STANDARD)
    e = p.edges
    p.Set(name='gjl-edges', edges=e)
    pickedRegions = (p.sets['gjl-edges'].edges, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))

    # 创建作业
    a.regenerate()
    mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1,
        multiprocessingMode=DEFAULT, numCpus=4, numDomains=4, numGPUs=0)

    # 可视图
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)