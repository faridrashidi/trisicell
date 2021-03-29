# Copyright Ian Smith 2002-2018
# Above only included as I am told that without it, others may claim copyright and stop me from editting my own copies of my own code! No idea how true this is.
# Beyond that, I am happy for anyone to make any edits they like to this code.
# For example, raising ValueError in pmf, cdf and sf functions should perhaps only be done when the distribution objects are created

# This code is basically translated from the VBA code available at https://iandjmsmith.wordpress.com

import math
import time

NonIntegralValuesAllowed_Others = False

OneOverSqrTwoPi = 0.39894228040143267793994605993438
cfSmall = 0.00000000000001
cfVSmall = 0.000000000000001
scalefactor = 1.1579208923731619542357098500869e77  # 2**256  ' used for rescaling calcs w/o impacting accuracy, to avoid over/underflow
scalefactor2 = 8.6361685550944446253863518628004e-78  # 2**-256
max_discrete = 9007199254740991.0
minLog1Value = -0.79149064

HalfMinusEulers_const = -0.07721566490153286060651209011
Onep25Minusln2Minuseulers_const = -0.020362845461478170023744211558177
coeffs0Minusp25 = 0.07246703342411321823620758332301

# For logfbit functions                      # Stieltjes' continued fraction
cf_0 = 1.0 / 12.0
cf_1 = 1.0 / 30.0
cf_2 = 53.0 / 210.0
cf_3 = 195.0 / 371.0
cf_4 = 22999.0 / 22737.0
cf_5 = 29944523.0 / 19733142.0
cf_6 = 109535241009.0 / 48264275462.0
cf_7 = 3.0099173832593981700731407342077
cf_8 = 4.026887192343901226168879531814
cf_9 = 5.0027680807540300516885024122767
cf_10 = 6.2839113708157821800726631549524
cf_11 = 7.4959191223840339297523547082674
cf_12 = 9.0406602343677266995311393604326
cf_13 = 10.489303654509482277188371304593
cf_14 = 12.297193610386205863989437140092
cf_15 = 13.982876953992430188259760651279
cf_16 = 16.053551416704935469715616365007
cf_17 = 17.976607399870277592569472307671
cf_18 = 20.309762027441653743805414720495
cf_19 = 22.470471639933132495517941571508
cf_20 = 25.065846548945972029163400322506
cf_21 = 27.464451825029133609175558982646
cf_22 = 30.321821231673047126882599306406
cf_23 = 32.958533929972987219994066451412
cf_24 = 36.077698931299242645153320900855
cf_25 = 38.952706682311555734544390410481
cf_26 = 42.333490043576957211381853948856
cf_27 = 45.446960850061621014424175737541
cf_28 = 49.089203129012597708164883350275
cf_29 = 52.441288751415337312569856046996
cf_30 = 56.344845345341843538420365947476
cf_31 = 59.935683907165858207852583492752
cf_32 = 64.100422755920354527906611892238
cf_33 = 67.930140788018221186367702745199
cf_34 = 72.355940555211701969680052963236
cf_35 = 76.424654626829689752585090422288
cf_36 = 81.111403237247965484814230985683
cf_37 = 85.419221276410972614585638717349
cf_38 = 90.366814723864108595513574581683
cf_39 = 94.913837100009887953076231291987
cf_40 = 100.12217846392919748899074683447

lfbArray = [
    cf_0,
    cf_1,
    cf_2,
    cf_3,
    cf_4,
    cf_5,
    cf_6,
    cf_7,
    cf_8,
    cf_9,
    cf_10,
    cf_11,
    cf_12,
    cf_13,
    cf_14,
    cf_15,
    cf_16,
    cf_17,
    cf_18,
    cf_19,
    cf_20,
    cf_21,
    cf_22,
    cf_23,
    cf_24,
    cf_25,
    cf_26,
    cf_27,
    cf_28,
    cf_29,
]
coeffs = [
    0.32246703342411321824,
    6.7352301053198095133e-02,
    2.0580808427784547879e-02,
    7.3855510286739852663e-03,
    2.8905103307415232858e-03,
    1.1927539117032609771e-03,
    5.0966952474304242234e-04,
    2.2315475845357937976e-04,
    9.9457512781808533715e-05,
    4.4926236738133141700e-05,
    2.0507212775670691553e-05,
    9.4394882752683959040e-06,
    4.3748667899074878042e-06,
    2.0392157538013662368e-06,
    9.5514121304074198329e-07,
    4.4924691987645660433e-07,
    2.1207184805554665869e-07,
    1.0043224823968099609e-07,
    4.7698101693639805658e-08,
    2.2711094608943164910e-08,
    1.0838659214896954091e-08,
    5.1834750419700466551e-09,
    2.4836745438024783172e-09,
    1.1921401405860912074e-09,
    5.7313672416788620133e-10,
    2.7595228851242331452e-10,
    1.3304764374244489481e-10,
    6.4229645638381000221e-11,
    3.1044247747322272762e-11,
    1.5021384080754142171e-11,
    7.2759744802390796625e-12,
    3.5277424765759150836e-12,
    1.7119917905596179086e-12,
    8.3153858414202848198e-13,
    4.0422005252894400655e-13,
    1.9664756310966164904e-13,
    9.5736303878385557638e-14,
    4.6640760264283742246e-14,
    2.2737369600659723206e-14,
    1.1091399470834522017e-14,
    5.4136591567253631315e-15,
    2.6438800178609949985e-15,
    1.2918959062789967293811764562316e-15,
    6.3159355041984485676779394847024e-16,
    1.421085482803160676983430580383e-14,
]
coeffs2 = [
    9.96703342411321823620758332301e-3,
    4.85230105319809513323333333333e-3,
    2.35164176111788121233425746862e-3,
    1.1355510286739852662730972914e-3,
    5.4676033074152328575298829848e-4,
    2.6269438789373716759012616902e-4,
    1.2601997117161385090708338501e-4,
    6.03943417869127130948e-5,
    2.89279988929196448257e-5,
    1.38537935563149598820e-5,
    6.63558635521614609862e-6,
    3.17947224962737026300e-6,
    1.52432377823166362836e-6,
    7.31319548444223379639e-7,
    3.51147479316783649952e-7,
    1.68754473874618369035e-7,
    8.11753732546888155551e-8,
    3.90847775936649142159e-8,
    1.88369052760822398681e-8,
    9.08717580313959348175e-9,
    4.38794008336117216467e-9,
    2.12078578473653627963e-9,
    1.02595225309999020577e-9,
    4.96752618206533915776e-10,
    2.40726205228207715756e-10,
    1.16751848407213311847e-10,
    5.66693373586358102002e-11,
    2.75272781658498271913e-11,
    1.33812334011666457419e-11,
    6.50929603319331702812e-12,
    3.16857905287746826547e-12,
    1.54339039998043529180e-12,
    7.52239805801019839357e-13,
    3.66855576849641617566e-13,
    1.79011840661361776213e-13,
    8.73989502840846835258e-14,
    4.26932273880725309600e-14,
    2.08656067727432658676e-14,
    1.02026669800712891581e-14,
    4.99113012967463749398e-15,
    2.44274876330334144841e-15,
    1.19604099925791156327e-15,
    5.85859784064943690566e-16,
    2.87087981566462948467e-16,
    1.40735520163755175636e-16,
]

logfbit0p5 = 0.054814121051917653896138702348386  # logfbit(0.5)
lfb_0 = 8.1061466795327258219670263594382e-02  # logfbit(0.0)
lfb_1 = 4.1340695955409294093822081407118e-02  # logfbit(1.0)
lfb_2 = 2.7677925684998339148789292746245e-02  # logfbit(2.0)
lfb_3 = 2.0790672103765093111522771767849e-02  # logfbit(3.0)
lfb_4 = 1.6644691189821192163194865373593e-02  # logfbit(4.0)
lfb_5 = 1.3876128823070747998745727023763e-02  # logfbit(5.0)

# For logfbit functions                   #Stirling's series for ln(Gamma(x)), A046968/A046969
lfbc1 = 1.0 / 12.0
lfbc2 = 1.0 / 30.0  # lfbc2 on are Sloane's ratio times 12
lfbc3 = 1.0 / 105.0
lfbc4 = 1.0 / 140.0
lfbc5 = 1.0 / 99.0
lfbc6 = 691.0 / 30030.0
lfbc7 = 1.0 / 13.0
lfbc8 = 0.35068485511628418514  # Chosen to make logfbit(6) & logfbit(7) correct
lfbc9 = 1.6769380337122674863  # Chosen to make logfbit(6) & logfbit(7) correct

# Global variables used by functions hypergeometric & CBNB0
ha1, hprob, hswap = 0.0, 0.0, False


def AlterForIntegralChecks_Others(value):
    # If non integral values are allowed then returns math.floor(value) otherwise only allows value if it equals math.floor(value)
    if NonIntegralValuesAllowed_Others:
        return math.floor(value)
    elif value != math.floor(value):
        raise ValueError
    else:
        return value


def Generalabminuscd(a, b, c, d):
    if (a == 0.0) or (b == 0.0) or (c == 0.0) or (d == 0.0):
        return a * b - c * d
    amr, ae = math.frexp(a)
    bmr, be = math.frexp(b)
    cmr, ce = math.frexp(c)
    dmr, de = math.frexp(d)
    am = int(amr * (2 ** 53))
    bm = int(bmr * (2 ** 53))
    cm = int(cmr * (2 ** 53))
    dm = int(dmr * (2 ** 53))
    abe = ae + be
    cde = ce + de
    abm = am * bm
    cdm = cm * dm
    scale = abe - cde
    if scale >= 0:
        if scale > 106:
            scale = 0
            cdm = 1
        # return math.ldexp(abm * (2**scale) - cdm, abe - scale - 106)
        return float(abm * (2 ** scale) - cdm) * (2.0 ** (abe - scale - 106))
    else:
        scale = -scale
        if scale > 106:
            scale = 0
            abm = 1
        # return math.ldexp(abm - cdm * (2**scale), cde - scale - 106)
        return float(abm - cdm * (2 ** scale)) * (2.0 ** (cde - scale - 106))


def logfbit2dif(x):
    # Calculation of logfbit2(x)-logfbit2(1+x).
    return 0.5 * (((x + 1.0) * (x + 2.0)) ** -2)


def logfbit2(x):
    # Second derivative of error part of Stirling's formula where log(x!) = log(sqrt(twopi))+(x+0.5)*log(x+1)-(x+1)+logfbit(x).
    if x >= 10000000000.0:
        return 2.0 * lfbc1 * ((x + 1.0) ** -3)
    elif x >= 7.0:
        x1 = x + 1.0
        x2 = 1.0 / (x1 * x1)
        x3 = x2 * (240.0 * lfbc8 - x2 * 306.0 * lfbc9)
        x3 = x2 * (132.0 * lfbc6 - x2 * (182.0 * lfbc7 - x3))
        x3 = x2 * (56.0 * lfbc4 - x2 * (90.0 * lfbc5 - x3))
        x3 = x2 * (12.0 * lfbc2 - x2 * (30.0 * lfbc3 - x3))
        return lfbc1 * (2.0 - x3) * x2 / x1
    elif x > -1.0:
        x1 = x
        x2 = 0.0
        while x1 < 7.0:
            x2 = x2 + logfbit2dif(x1)
            x1 = x1 + 1.0
        return x2 + logfbit2(x1)
    else:
        raise ValueError


def logfbit4dif(x):
    # Calculation of logfbit4(x)-logfbit4(1+x).
    return (10.0 * x * (x + 3.0) + 23.0) * (((x + 1.0) * (x + 2.0)) ** -4)


def logfbit4(x):
    # Fourth derivative of error part of Stirling's formula where log(x!) = log(sqrt(twopi))+(x+0.5)*log(x+1)-(x+1)+logfbit(x).
    if x >= 10000000000.0:
        return -0.5 * ((x + 1.0) ** -4)
    elif x >= 7.0:
        x1 = x + 1.0
        x2 = 1.0 / (x1 * x1)
        x3 = x2 * (73440.0 * lfbc8 - x2 * 116280.0 * lfbc9)
        x3 = x2 * (24024.0 * lfbc6 - x2 * (43680.0 * lfbc7 - x3))
        x3 = x2 * (5040.0 * lfbc4 - x2 * (11880.0 * lfbc5 - x3))
        x3 = x2 * (360.0 * lfbc2 - x2 * (1680.0 * lfbc3 - x3))
        return lfbc1 * (24.0 - x3) * x2 * x2 / x1
    elif x > -1.0:
        x1 = x
        x2 = 0.0
        while x1 < 7.0:
            x2 = x2 + logfbit4dif(x1)
            x1 = x1 + 1.0

        return x2 + logfbit4(x1)
    else:
        raise ValueError


def log0(x):
    # Accurate and quicker calculation of log(1+x), particularly for small x. Code from Wolfgang Ehrhardt.
    if x > 4.0:
        return math.log(1.0 + x)
    else:
        y = 1.0 + x
        if y == 1.0:
            return x
        else:
            return math.log(y) + (x - (y - 1.0)) / y


def logcf(x, i, d):
    # Continued fraction for calculation of 1/i + x/(i+d) + x*x/(i+2*d) + x*x*x/(i+3d) + ...
    c1 = 2.0 * d
    c2 = i + d
    c4 = c2 + d
    a1 = c2
    b1 = i * (c2 - i * x)
    b2 = d * d * x
    a2 = c4 * c2 - b2
    b2 = c4 * b1 - i * b2

    while abs(a2 * b1 - a1 * b2) > abs(cfVSmall * b1 * a2):

        c3 = c2 * c2 * x
        c2 = c2 + d
        c4 = c4 + d
        a1 = c4 * a2 - c3 * a1
        b1 = c4 * b2 - c3 * b1

        c3 = c1 * c1 * x
        c1 = c1 + d
        c4 = c4 + d
        a2 = c4 * a1 - c3 * a2
        b2 = c4 * b1 - c3 * b2
        if b2 > scalefactor:
            a1 = a1 * scalefactor2
            b1 = b1 * scalefactor2
            a2 = a2 * scalefactor2
            b2 = b2 * scalefactor2
        elif b2 < scalefactor2:
            a1 = a1 * scalefactor
            b1 = b1 * scalefactor
            a2 = a2 * scalefactor
            b2 = b2 * scalefactor
    return a2 / b2


def log1(x):
    # Accurate calculation of log(1+x)-x, particularly for small x.
    if abs(x) < 0.01:
        term = x / (2.0 + x)
        y = term * term
        return term * (
            (((2.0 / 9.0 * y + 2.0 / 7.0) * y + 0.4) * y + 2.0 / 3.0) * y - x
        )
    elif x < minLog1Value or x > 1.0:
        return math.log(1.0 + x) - x
    else:
        term = x / (2.0 + x)
        y = term * term
        return term * (2.0 * y * logcf(y, 3.0, 2.0) - x)


def logfbitdif(x):
    # Calculation of logfbit(x)-logfbit(1+x). x must be > -1.
    if x < -0.65:
        return (x + 1.5) * log0(1.0 / (x + 1.0)) - 1.0
    else:
        y2 = (2.0 * x + 3.0) ** -2
        return y2 * logcf(y2, 3.0, 2.0)


def logfbita(x):
    # Error part of Stirling's formula where math.log(x!) = math.log(math.sqrtt(twopi))+(x+0.5)*math.log(x+1)-(x+1)+logfbita(x).
    # Are we ever concerned about the relative error involved in this function? I don't think so.
    if x >= 100000000.0:
        return lfbc1 / (x + 1.0)
    elif x >= 6.0:  # Abramowitz & Stegun's series 6.1.41
        x1 = x + 1.0
        x2 = 1.0 / (x1 * x1)
        x3 = x2 * (lfbc6 - x2 * (lfbc7 - x2 * (lfbc8 - x2 * lfbc9)))
        x3 = x2 * (lfbc4 - x2 * (lfbc5 - x3))
        x3 = x2 * (lfbc2 - x2 * (lfbc3 - x3))
        return lfbc1 * (1.0 - x3) / x1
    elif x == 0.0:
        return lfb_0
    elif x == 1.0:
        return lfb_1
    elif x == 2.0:
        return lfb_2
    elif x == 3.0:
        return lfb_3
    elif x == 4.0:
        return lfb_4
    elif x == 5.0:
        return lfb_5
    elif x > -1.0:
        x1 = x
        x2 = 0.0
        while x1 < 6.0:
            x2 = x2 + logfbitdif(x1)
            x1 = x1 + 1.0
        return x2 + logfbit(x1)
    else:
        raise ValueError


def logfbit(x):
    # Calculates log of x factorial - log(sqrt(2*pi)) +(x+1) -(x+0.5)*log(x+1)
    # using the error part of Stirling's formula (see Abramowitz & Stegun# s series 6.1.41)
    # and Stieltjes' continued fraction for the gamma function.
    # For x < 1.5, uses expansion of log(x!) and log((x+1)!) from Abramowitz & Stegun# s series 6.1.33
    # We are primarily concerned about the absolute error in this function.
    # Due to cancellation errors in calculating 1+x as x tends to -1, the function loses accuracy and should not be used!
    if x >= 6.0:
        x1 = x + 1.0
        if x >= 1000.0:
            if x >= 100000000.0:
                x3 = 0.0
            else:
                x2 = 1.0 / (x1 * x1)
                x3 = x2 * (lfbc2 - x2 * lfbc3)

        else:
            x2 = 1.0 / (x1 * x1)
            if x >= 40.0:
                x3 = 0.0
            elif x >= 15.0:
                x3 = x2 * (lfbc6 - x2 * lfbc7)
            else:
                x3 = x2 * (lfbc6 - x2 * (lfbc7 - x2 * (lfbc8 - x2 * lfbc9)))

            x3 = x2 * (lfbc4 - x2 * (lfbc5 - x3))
            x3 = x2 * (lfbc2 - x2 * (lfbc3 - x3))

        return lfbc1 * (1.0 - x3) / x1
        # return (1.0 - x3) / (12.0 * x1)
    elif x == 0.0:
        return lfb_0
    elif x == 1.0:
        return lfb_1
    elif x == 2.0:
        return lfb_2
    elif x == 3.0:
        return lfb_3
    elif x == 4.0:
        return lfb_4
    elif x == 5.0:
        return lfb_5
    elif x > 1.5:
        x1 = x + 1.0
        if x >= 2.5:
            # x2 = 0.25 * ((math.sqrt(x1 * x1 + 81.0) - x1) + 81.0 / (x1 + math.sqrt(x1 * x1 + 90.25)))
            x2 = 40.5 / (x1 + math.sqrt(x1 * x1 + 81.0))
        else:
            # x2 = 0.25 * ((math.sqrt(x1 * x1 + 225.0) - x1) + 225.0 / (x1 + math.sqrt(x1 * x1 + 240.25)))
            x2 = 112.5 / (x1 + math.sqrt(x1 * x1 + 225.0))
            x2 = cf_27 / (x1 + cf_28 / (x1 + cf_29 / (x1 + x2)))
            x2 = cf_24 / (x1 + cf_25 / (x1 + cf_26 / (x1 + x2)))
            x2 = cf_21 / (x1 + cf_22 / (x1 + cf_23 / (x1 + x2)))
            x2 = cf_18 / (x1 + cf_19 / (x1 + cf_20 / (x1 + x2)))

        x2 = cf_15 / (x1 + cf_16 / (x1 + cf_17 / (x1 + x2)))
        x2 = cf_12 / (x1 + cf_13 / (x1 + cf_14 / (x1 + x2)))
        x2 = cf_9 / (x1 + cf_10 / (x1 + cf_11 / (x1 + x2)))
        x2 = cf_6 / (x1 + cf_7 / (x1 + cf_8 / (x1 + x2)))
        x2 = cf_3 / (x1 + cf_4 / (x1 + cf_5 / (x1 + x2)))
        # return cf_0 / (x1 + cf_1 / (x1 + cf_2 / (x1 + x2)))
        return 1.0 / (12.0 * (x1 + cf_1 / (x1 + cf_2 / (x1 + x2))))
    # elif (x == 1.5):
    #    return 3.316287351993628748511050974106e-02   #  3.316287351993628748511050974106e-02
    elif x == 0.5:
        return logfbit0p5  #  5.481412105191765389613870234839e-02
    elif x == -0.5:
        return 0.15342640972002734529138393927091  #  0.15342640972002734529138393927091
    elif x >= -0.65:
        if x <= 0.0:
            i = len(coeffs) - 1
            lgam = coeffs[i] * logcf(-x / 2.0, i + 2.0, 1.0)
            # print(i, coeffs[i],lgam)
            for i in range(i - 1, 0, -1):
                lgam = coeffs[i] - x * lgam
                # print(i,lgam)
            return (
                (
                    coeffs0Minusp25
                    - (x * (1.0 / 3.0 - (x + 1.5) * logcf(-x, 4.0, 1.0)) + lgam) * x
                )
                * x
                + HalfMinusEulers_const
            ) * x + lfb_0
        elif x <= 1.56:
            x = x - 1.0
            i = len(coeffs2) + 2
            lgam = ((x + 2.5) * logcf(-x / 2.0, i, 1.0) - (2.0 / (i - 1.0))) * (
                2.0 ** -i
            ) + (3.0 ** -i) * logcf(-x / 3.0, i, 1.0)
            # print(i,coeffs2[i-3],lgam)
            for i in range(i - 3, -1, -1):
                lgam = coeffs2[i] - x * lgam
                # print(i,lgam)
            return (x * lgam + Onep25Minusln2Minuseulers_const) * x + lfb_1
        elif x <= 2.5:
            x = x - 2.0
            i = len(coeffs) - 1
            lgam = coeffs[i] * logcf(-x / 2.0, i + 2.0, 1.0)
            # print(i, coeffs[i],lgam)
            for i in range(i - 1, 0, -1):
                lgam = coeffs[i] - x * lgam
                # print(lgam)
            return (
                (
                    x * x * (coeffs0Minus1Third - x * lgam)
                    - (x + 2.5) * log1(x / 3.0)
                    + log1(x / 2.0)
                )
                + FiveOver3Minusln3Minuseulers_const * x
            ) + lfb_2
        else:
            x = x - 3.0
            i = len(coeffs) - 1
            lgam = coeffs[i] * logcf(-x / 2.0, i + 2.0, 1.0)
            # print(i, coeffs[i],lgam)
            for i in range(i - 1, 0, -1):
                lgam = coeffs[i] - x * lgam
                # print(lgam)
            return (
                (
                    x * x * (coeffs0Minusp25 - x * lgam)
                    - (x + 3.5) * log1(x / 4.0)
                    + log1(x / 2.0)
                    + log1(x / 3.0)
                )
                + Forty7Over48Minusln4Minuseulers_const * x
            ) + lfb_3

    elif x > -1.0:
        return logfbitdif(x) + logfbit(x + 1.0)
    else:
        raise ValueError


def lfbaccdif1(a, b):
    # Calculates logfbit(b)-logfbit(a+b) accurately for a > 0 & b >= 0. Reasonably accurate for a >=0 & b < 0.
    if a < 0.0:
        return -lfbaccdif1(-a, b + a)
    elif b >= 8.0:
        y1 = b + 1.0
        y2 = y1 ** -2
        x1 = a + b + 1.0
        x2 = x1 ** -2
        x3 = x2 * lfbc9
        y3 = y2 * lfbc9
        acc = x2 * (a * (x1 + y1) * y3)
        x3 = x2 * (lfbc8 - x3)
        y3 = y2 * (lfbc8 - y3)
        acc = x2 * (a * (x1 + y1) * y3 - acc)
        x3 = x2 * (lfbc7 - x3)
        y3 = y2 * (lfbc7 - y3)
        acc = x2 * (a * (x1 + y1) * y3 - acc)
        x3 = x2 * (lfbc6 - x3)
        y3 = y2 * (lfbc6 - y3)
        acc = x2 * (a * (x1 + y1) * y3 - acc)
        x3 = x2 * (lfbc5 - x3)
        y3 = y2 * (lfbc5 - y3)
        acc = x2 * (a * (x1 + y1) * y3 - acc)
        x3 = x2 * (lfbc4 - x3)
        y3 = y2 * (lfbc4 - y3)
        acc = x2 * (a * (x1 + y1) * y3 - acc)
        x3 = x2 * (lfbc3 - x3)
        y3 = y2 * (lfbc3 - y3)
        acc = x2 * (a * (x1 + y1) * y3 - acc)
        x3 = x2 * (lfbc2 - x3)
        y3 = y2 * (lfbc2 - y3)
        acc = x2 * (a * (x1 + y1) * y3 - acc)
        # return lfbc1 * (a * (1.0 - y3) - y1 * acc) / (x1 * y1)
        return (a * (1.0 - y3) - y1 * acc) / (12.0 * x1 * y1)
    elif b >= 1.7:
        y1 = b + 1.0
        x1 = a + b + 1.0
        if b >= 3.0:
            Start = 17
        else:
            Start = 29

        s1 = (0.5 * (Start + 1.0)) ** 2
        s2 = (0.5 * (Start + 1.5)) ** 2
        ty = y1 * math.sqrt(1.0 + s1 * (y1 ** -2))
        tx = x1 * math.sqrt(1.0 + s1 * (x1 ** -2))
        y2 = ty - y1
        x2 = tx - x1
        acc = a * (1.0 - (2.0 * y1 + a) / (tx + ty))
        # Seems to work better without the next 2 lines. - Not with modification to s2
        ty = y1 * math.sqrt(1.0 + s2 * (y1 ** -2))
        tx = x1 * math.sqrt(1.0 + s2 * (x1 ** -2))
        acc = 0.25 * (
            acc + s1 / ((y1 + ty) * (x1 + tx)) * a * (1.0 + (2.0 * y1 + a) / (tx + ty))
        )
        y2 = 0.25 * (y2 + s1 / (y1 + ty))
        x2 = 0.25 * (x2 + s1 / (x1 + tx))
        for i in range(Start, 0, -1):
            acc = lfbArray[i] * (a - acc) / ((x1 + x2) * (y1 + y2))
            y2 = lfbArray[i] / (y1 + y2)
            x2 = lfbArray[i] / (x1 + x2)
        return cf_0 * (a - acc) / ((x1 + x2) * (y1 + y2))
        # return (a - acc) / (12.0 * (x1 + x2) * (y1 + y2))
    elif b > -1.0:
        if b < -0.66:
            if a > 1.0:
                return logfbitdif(b) + lfbaccdif1(a - 1.0, b + 1.0)
            elif a == 1.0:
                return logfbitdif(b)
            else:
                s2 = a * log0(1.0 / (b + 1.0 + a))
                s1 = logfbitdif(b + a)
                if s1 > s2:
                    s1 = (b + 1.5) * log0(a / ((b + 1.0) * (b + 2.0 + a))) - s2
                else:
                    s2 = s1
                    s1 = logfbitdif(b) - s1

                if s1 > 0.1 * s2:
                    return s1 + lfbaccdif1(a, b + 1.0)

        if b + a > 2:
            s1 = lfbaccdif1(b + a - 1.75, 1.75)
            a = 1.75 - b
        else:
            s1 = 0.0

        y1 = b - 1.0
        x1 = y1 + a
        i = len(coeffs2) + 2
        scale2 = 2.0 ** -i
        scale3 = 3.0 ** -i
        # y2 = ((y1 + 2.5) * logcf(-y1 / 2.0, i, 1.0) - (2.0 / (i - 1.0))) * scale2 + (scale3 * logcf(-y1 / 3.0, i, 1.0) + scale2 * scale2 * logcf(-y1 / 4.0, i, 1.0))
        # x2 = ((x1 + 2.5) * logcf(-x1 / 2.0, i, 1.0) - (2.0 / (i - 1.0))) * scale2 + (scale3 * logcf(-x1 / 3.0, i, 1.0) + scale2 * scale2 * logcf(-x1 / 4.0, i, 1.0))
        y2 = (
            (y1 + 2.5) * logcf(-y1 / 2.0, i, 1.0) - (2.0 / (i - 1.0))
        ) * scale2 + scale3 * logcf(-y1 / 3.0, i, 1.0)
        x2 = (
            (x1 + 2.5) * logcf(-x1 / 2.0, i, 1.0) - (2.0 / (i - 1.0))
        ) * scale2 + scale3 * logcf(-x1 / 3.0, i, 1.0)
        if a > 0.000006:
            acc = (
                y2 - x2
            )  # This calculation is not accurate enough for b < 0 and a small - hence if b < 0 code above and derivative code below for small a
        else:
            y3 = -(y1 + a / 2.0) / 2.0
            x3 = -(y1 + a / 2.0) / 3.0
            acc = -a * (
                scale2
                * (
                    logcf(y3, i, 1.0)
                    + (y3 - 1.25) * (1.0 / (1.0 - y3) - i * logcf(y3, i + 1.0, 1.0))
                )
                - scale3 / 3.0 * ((1.0 / (1.0 - x3) - i * logcf(x3, i + 1.0, 1.0)))
            )

        for i in range(i - 3, -1, -1):
            acc = a * y2 - x1 * acc
            y2 = coeffs2[i] - y1 * y2
            x2 = coeffs2[i] - x1 * x2
        return s1 + (
            y1 * y1 * acc - a * (x2 * (x1 + y1) + Onep25Minusln2Minuseulers_const)
        )
    else:
        raise ValueError


def hypergeometricTerm(ai, aji, aki, amkji):
    # Probability that hypergeometric variate from a population with total type Is of aki+ai, total type IIs of amkji+aji, has ai type Is and aji type IIs selected.
    # Parameterised this way for use with the Beta Negative Binomial distribution.
    ak = aki + ai
    amk = amkji + aji
    aj = aji + ai
    am = amk + ak
    amj = amkji + aki
    if am > max_discrete:
        raise ValueError
    if (ai == 0.0) and ((aji <= 0.0) or (aki <= 0.0) or (amj < 0.0) or (amk < 0.0)):
        return 1.0
    elif (ai > 0.0) and (min(aki, aji) == 0.0) and (max(amj, amk) == 0.0):
        return 1.0
    elif (ai >= 0.0) and (amkji > -1.0) and (aki > -1.0) and (aji >= 0.0):
        c1 = logfbit(amkji) + logfbit(aki) + logfbit(aji) + logfbit(am) + logfbit(ai)
        c1 = logfbit(amk) + logfbit(ak) + logfbit(aj) + logfbit(amj) - c1
        # c1 = lfbaccdif1(ak, amk) - lfbaccdif1(ai, aki) - lfbaccdif1(ai, aji) - lfbaccdif1(aki, amkji) - logfbit(ai)
        ai1 = ai + 1.0
        aj1 = aj + 1.0
        ak1 = ak + 1.0
        am1 = am + 1.0
        aki1 = aki + 1.0
        aji1 = aji + 1.0
        amk1 = amk + 1.0
        amj1 = amj + 1.0
        amkji1 = amkji + 1.0
        cjkmi = Generalabminuscd(aji, aki, ai, amkji)
        c5 = (cjkmi - ai) / (amkji1 * am1)
        if c5 < minLog1Value:
            c3 = amkji * (math.log((amj1 * amk1) / (amkji1 * am1)) - c5) - c5
        else:
            c3 = amkji * log1(c5) - c5

        c5 = (-cjkmi - aji) / (aki1 * am1)
        if c5 < minLog1Value:
            c4 = aki * (math.log((ak1 * amj1) / (aki1 * am1)) - c5) - c5
        else:
            c4 = aki * log1(c5) - c5

        c3 = c3 + c4
        c5 = (-cjkmi - aki) / (aji1 * am1)
        if c5 < minLog1Value:
            c4 = aji * (math.log((aj1 * amk1) / (aji1 * am1)) - c5) - c5
        else:
            c4 = aji * log1(c5) - c5

        c3 = c3 + c4
        c5 = (cjkmi - amkji) / (ai1 * am1)
        if c5 < minLog1Value:
            c4 = ai * (math.log((aj1 * ak1) / (ai1 * am1)) - c5) - c5
        else:
            c4 = ai * log1(c5) - c5

        c3 = c3 + c4

        return (
            math.exp((c1 + 1.0 / am1) + c3)
            * math.sqrt(
                (amk1 * ak1) * (aj1 * amj1) / ((amkji1 * aki1 * aji1) * (am1 * ai1))
            )
            * OneOverSqrTwoPi
        )
    else:
        return 0.0


def hypergeometric(ai, aji, aki, amkji, comp):  # , ha1, hprob, hswap):
    # Probability that hypergeometric variate from a population with total type Is of aki+ai, total type IIs of amkji+aji, has up to ai type Is selected in a sample of size aji+ai.
    # Parameterised this way for use with the Beta Negative Binomial distribution.
    global ha1, hprob, hswap
    if (amkji > -1.0) and (amkji < 0.0):
        ip1 = -amkji
        mkji = ip1 - 1.0
        allintegral = False
    else:
        ip1 = amkji + 1.0
        mkji = amkji
        allintegral = (
            ai == math.floor(ai)
            and aji == math.floor(aji)
            and aki == math.floor(aki)
            and mkji == math.floor(mkji)
        )

    if allintegral:
        swapped = (ai + 0.5) * (mkji + 0.5) >= (aki - 0.5) * (aji - 0.5)
    elif ai < 100.0 and ai == math.floor(ai) or mkji < 0.0:
        if comp:
            swapped = (ai + 0.5) * (mkji + 0.5) >= aki * aji
        else:
            swapped = (ai + 0.5) * (mkji + 0.5) >= aki * aji + 1000.0

    elif ai < 1.0:
        swapped = (ai + 0.5) * (mkji + 0.5) >= aki * aji
    elif aji < 1.0 or aki < 1.0 or (ai < 1.0 and ai > 0.0):
        swapped = False
    else:
        swapped = (ai + 0.5) * (mkji + 0.5) >= (aki - 0.5) * (aji - 0.5)

    if not swapped:
        i = ai
        ji = aji
        ki = aki
    else:
        i = aji - 1.0
        ji = ai + 1.0
        ki = ip1
        ip1 = aki
        mkji = aki - 1.0

    c2 = ji + i
    c4 = mkji + ki + c2
    if c4 > max_discrete:
        raise ValueError

    if (
        (i >= 0.0)
        and (ji <= 0.0)
        or (ki <= 0.0)
        or (ip1 + ki <= 0.0)
        or (ip1 + ji <= 0.0)
    ):
        exact = True
        if i >= 0.0:
            prob = 1.0
        else:
            prob = 0.0

    elif (ip1 > 0.0) and (ip1 < 1.0):
        exact = False
        prob = (
            hypergeometricTerm(i, ji, ki, ip1)
            * (ip1 * (c4 + 1.0))
            / ((ki + ip1) * (ji + ip1))
        )
    else:
        exact = (
            (i == 0.0)
            and ((ji <= 0.0) or (ki <= 0.0) or (mkji + ki < 0.0) or (mkji + ji < 0.0))
        ) or ((i > 0.0) and (min(ki, ji) == 0.0) and (max(mkji + ki, mkji + ji) == 0.0))
        prob = hypergeometricTerm(i, ji, ki, mkji)

    hprob = prob
    hswap = swapped
    ha1 = 0.0

    if (exact) or (prob == 0.0):
        if swapped == comp:
            return prob
        else:
            return 1.0 - prob

    a1 = 0.0
    sumAlways = 0.0
    sumFactor = 10.0

    if i < mkji:
        must_do_cf = i != math.floor(i)
        maxSums = math.floor(i)
    else:
        must_do_cf = mkji != math.floor(mkji)
        maxSums = math.floor(max(mkji, 0.0))

    if must_do_cf:
        sumAlways = 0.0
        sumFactor = 5.0
    else:
        sumAlways = 20.0
        sumFactor = 10.0

    if maxSums > sumAlways or must_do_cf:
        numb = math.floor(
            sumFactor
            / c4
            * math.exp(math.log((ki + i) * (ji + i) * (ip1 + ji) * (ip1 + ki)) / 3.0)
        )
        numb = math.floor(i - (ki + i) * (ji + i) / c4 + numb)
        if numb < 0.0:
            numb = 0.0
        elif numb > maxSums:
            numb = maxSums
    else:
        numb = maxSums

    if 2.0 * numb <= maxSums or must_do_cf:
        b1 = 1.0
        c1 = 0.0
        c2 = i - numb
        c3 = mkji - numb
        s = c3
        a2 = c2
        c3 = c3 - 1.0
        b2 = Generalabminuscd(ki + numb + 1.0, ji + numb + 1.0, c2 - 1.0, c3)
        bn = b2
        bnAdd = c3 + c4 + c2 - 2.0
        while b2 > 0.0 and (abs(a2 * b1 - a1 * b2) > abs(cfVSmall * b1 * a2)):
            c1 = c1 + 1.0
            c2 = c2 - 1.0
            an = (c1 * c2) * (c3 * c4)
            c3 = c3 - 1.0
            c4 = c4 - 1.0
            bn = bn + bnAdd
            bnAdd = bnAdd - 4.0
            a1 = bn * a2 + an * a1
            b1 = bn * b2 + an * b1
            if b1 > scalefactor:
                a1 = a1 * scalefactor2
                b1 = b1 * scalefactor2
                a2 = a2 * scalefactor2
                b2 = b2 * scalefactor2
            c1 = c1 + 1.0
            c2 = c2 - 1.0
            an = (c1 * c2) * (c3 * c4)
            c3 = c3 - 1.0
            c4 = c4 - 1.0
            bn = bn + bnAdd
            bnAdd = bnAdd - 4.0
            a2 = bn * a1 + an * a2
            b2 = bn * b1 + an * b2
            if b2 > scalefactor:
                a1 = a1 * scalefactor2
                b1 = b1 * scalefactor2
                a2 = a2 * scalefactor2
                b2 = b2 * scalefactor2

        if b1 < 0.0 or b2 < 0.0:
            raise ValueError  # Actually, something unexpected has happened here - although it will be to do with the values input.
        else:
            a1 = a2 / b2 * s

    else:
        numb = maxSums

    c1 = i - numb + 1.0
    c2 = mkji - numb + 1.0
    c3 = ki + numb
    c4 = ji + numb
    while numb > 0.0:
        a1 = (1.0 + a1) * ((c1 * c2) / (c3 * c4))
        c1 = c1 + 1.0
        c2 = c2 + 1.0
        c3 = c3 - 1.0
        c4 = c4 - 1.0
        numb = numb - 1.0

    ha1 = a1
    a1 = (1.0 + a1) * prob
    if swapped == comp:
        return a1
    else:
        if a1 > 0.99:
            raise ValueError  # Similarly, if we end up here then something unexpected has happened - although it will be to do with the values input.
        else:
            return 1.0 - a1


def PBB(i, ssmi, beta_shape1, beta_shape2):
    global hTerm
    hTerm = hypergeometricTerm(i, ssmi, beta_shape1, beta_shape2)
    return (
        (beta_shape1 / (i + beta_shape1))
        * (beta_shape2 / (beta_shape1 + beta_shape2))
        * ((i + ssmi + beta_shape1 + beta_shape2) / (ssmi + beta_shape2))
        * hTerm
    )


def pmf_BetaNegativeBinomial(i, r, beta_shape1, beta_shape2):
    i = AlterForIntegralChecks_Others(i)
    if r <= 0.0 or beta_shape1 <= 0.0 or beta_shape2 <= 0.0:
        raise ValueError
    elif i < 0:
        return 0.0
    else:
        return (
            (beta_shape2 / (beta_shape1 + beta_shape2))
            * (r / (beta_shape1 + r))
            * beta_shape1
            * (i + beta_shape1 + r + beta_shape2)
            / ((i + r) * (i + beta_shape2))
            * hypergeometricTerm(i, r, beta_shape2, beta_shape1)
        )


def CBNB0(i, r, beta_shape1, beta_shape2, toBeAdded):
    global ha1, hprob, hswap
    if r < 2.0 or beta_shape2 < 2.0:
        # Assumption here that i is integral or greater than 4.
        mrb2 = max(r, beta_shape2)
        other = min(r, beta_shape2)
        cbnb0 = PBB(i, other, mrb2, beta_shape1)
        if i == 0.0:
            return cbnb0
        cbnb0 = cbnb0 * (
            1.0 + i * (other + beta_shape1) / (((i - 1.0) + mrb2) * (other + 1.0))
        )
        if i == 1.0:
            return cbnb0
        i = i - 2.0
        other = other + 2.0
        temp = PBB(i, mrb2, other, beta_shape1)
        if i == 0.0:
            return cbnb0 + temp
        cbnb0 = cbnb0 + temp * (
            1.0 + i * (mrb2 + beta_shape1) / (((i - 1.0) + other) * (mrb2 + 1.0))
        )
        if i == 1.0:
            return cbnb0
        i = i - 2.0
        mrb2 = mrb2 + 2.0
        return cbnb0 + CBNB0(i, mrb2, beta_shape1, other, cbnb0)
    elif beta_shape1 < 1.0:
        mrb2 = max(r, beta_shape2)
        other = min(r, beta_shape2)
        cbnb0 = hypergeometric(
            i, mrb2 - 1.0, other, beta_shape1, False
        )  # , ha1, hprob, hswap)
        if hswap:
            temp = PBB(mrb2 - 1.0, beta_shape1, i + 1.0, other)
            if (toBeAdded + (cbnb0 - temp)) < 0.01 * (toBeAdded + (cbnb0 + temp)):
                return CBNB2(i, mrb2, beta_shape1, other)
            else:
                return cbnb0 - temp

        elif ha1 < -0.9 * beta_shape1 / (beta_shape1 + other):
            raise ValueError
        else:
            return hprob * (beta_shape1 / (beta_shape1 + other) + ha1)
    else:
        return hypergeometric(
            i, r, beta_shape2, beta_shape1 - 1.0, False
        )  # , ha1, hprob, hswap)


def CBNB2(i, r, beta_shape1, beta_shape2):
    ss = min(r, beta_shape2)
    bs2 = max(r, beta_shape2)
    r = ss
    beta_shape2 = bs2
    d1 = (i + 0.5) * (beta_shape1 + 0.5) - (bs2 - 1.5) * (ss - 0.5)
    if d1 < 0.0:
        return CBNB0(i, ss, beta_shape1, bs2, 0.0)

    d1 = math.floor(d1 / (bs2 + beta_shape1 - 1.0)) + 10.0
    if ss + d1 > bs2:
        d1 = math.floor(bs2 - ss)
    ss = ss + d1
    j = i - d1
    d2 = (j + 0.5) * (beta_shape1 + 0.5) - (bs2 - 1.5) * (ss - 0.5)
    if d2 < 0.0:
        d2 = 10.0
    else:
        temp = bs2 + ss + 2.0 * beta_shape1 - 1.0
        d2 = math.floor((math.sqrt(temp * temp + 4.0 * d2) - temp) / 2.0) + 10.0

    if 2.0 * d2 > i:
        d2 = math.floor(i / 2.0)

    pbbval = PBB(i, r, beta_shape2, beta_shape1)
    ss = ss + d2
    bs2 = bs2 + d2
    j = j - 2.0 * d2
    cbnb2 = CBNB0(j, ss, beta_shape1, bs2, 0.0)
    temp = 1.0
    d_count = d2 - 2.0
    # print(d1,d2)
    j = j + 1.0
    while d_count >= 0.0:
        j = j + 1.0
        bs2 = beta_shape2 + d_count
        d_count = d_count - 1.0
        temp = 1.0 + (j * (bs2 + beta_shape1) / ((j + ss - 1.0) * (bs2 + 1.0))) * temp

    j = i - d2 - d1
    temp = (ss * (j + bs2)) / (bs2 * (j + ss)) * temp
    d_count = d1 + d2 - 1.0
    while d_count >= 0:
        j = j + 1.0
        ss = r + d_count
        d_count = d_count - 1.0
        temp = 1.0 + (j * (ss + beta_shape1) / ((j + bs2 - 1.0) * (ss + 1.0))) * temp

    return cbnb2 + temp * pbbval


def ccBNB5(ilim, rr, a, bb):
    if rr > bb:
        r = rr
        b = bb
    else:
        r = bb
        b = rr

    ccbnb5 = (
        (a + 0.5) * log0(b * r / ((a + 1.0) * (b + a + r + 1.0)))
        - r * log0(b / (r + a + 1.0))
        - b * log0(r / (a + b + 1.0))
    )
    if r <= 0.001:
        temp = a + (b + r) * 0.5
        ccbnb5 = ccbnb5 - b * r * (
            logfbit2(temp) + (b ** 2 + r ** 2) * logfbit4(temp) / 24.0
        )
    else:
        ccbnb5 = ccbnb5 + (lfbaccdif1(b, r + a) - lfbaccdif1(b, a))

    temp = 0.0
    if ilim > 0.0:
        i = ilim
        while i > 1.0:
            i = i - 1.0
            temp = (1.0 + temp) * (i + r) * (i + b) / ((i + r + a + b) * (i + 1.0))

        temp = (1.0 + temp) * math.exp(ccbnb5) * a

    return (r * b * (1.0 - temp) - math.expm1(ccbnb5) * a * (r + a + b)) / (
        (r + a) * (a + b)
    )


def cdf_BetaNegativeBinomial(i, r, beta_shape1, beta_shape2):
    i = math.floor(i)
    if r <= 0.0 or beta_shape1 <= 0.0 or beta_shape2 <= 0.0:
        raise ValueError
    elif i < 0:
        return 0.0
    else:
        return CBNB0(i, r, beta_shape1, beta_shape2, 0.0)


def sf_BetaNegativeBinomial(i, r, beta_shape1, beta_shape2):
    i = math.floor(i)
    mrb2 = max(r, beta_shape2)
    other = min(r, beta_shape2)
    if other <= 0.0 or beta_shape1 <= 0.0 or beta_shape2 <= 0.0:
        raise ValueError
    elif i < 0.0:
        return 1.0
    elif (i == 0.0) or (
        (i < 1000000.0)
        and (other < 0.001)
        and (beta_shape1 > 50.0 * other)
        and (10.0 * i * beta_shape1 < mrb2)
    ):
        return ccBNB5(i, mrb2, beta_shape1, other)
    elif (
        mrb2 >= 100.0
        or other > 20.0
        or (
            mrb2 >= 5.0
            and (other - 0.5) * (mrb2 - 0.5) > (i + 0.5) * (beta_shape1 + 0.5)
        )
    ):
        return CBNB0(mrb2 - 1.0, i + 1.0, other, beta_shape1, 0.0)
    else:
        ccbnb = 0.0
        temp = 0.0
        i = i + 1.0
        if other >= 1.0:
            mrb2 = mrb2 - 1.0
            other = other - 1.0
            temp = hypergeometricTerm(i, mrb2, other, beta_shape1)
            ccbnb = temp
            while (other >= 1.0) and (temp > 1e-16 * ccbnb):
                i = i + 1.0
                beta_shape1 = beta_shape1 + 1.0
                temp = temp * (mrb2 * other) / (i * beta_shape1)
                mrb2 = mrb2 - 1.0
                other = other - 1.0
                ccbnb = ccbnb + temp

            if other >= 1.0:
                return ccbnb
            i = i + 1.0
            beta_shape1 = beta_shape1 + 1.0

        if mrb2 >= 1.0:
            mxib1 = max(i, beta_shape1)
            mnib1 = min(i, beta_shape1)
            if temp == 0.0:
                mrb2 = mrb2 - 1.0
                temp = PBB(mnib1, mrb2, other, mxib1)
            else:  # temp is hypergeometricTerm(mnib1-1, mrb2, other, mxib1-1)
                temp = temp * other * mrb2
                mrb2 = mrb2 - 1.0
                temp = temp / (mnib1 * (mrb2 + mxib1))

            ccbnb = ccbnb + temp
            while (mrb2 >= 1.0) and (temp > 1e-16 * ccbnb):
                temp = temp * mrb2 * (mnib1 + other)
                mnib1 = mnib1 + 1.0
                if mnib1 > mxib1:
                    mxib1, mnib1 = mnib1, mxib1

                # Block below not required if hypergeometric block included above and therefore other guaranteed < 1 <= mrb2
                # if mrb2 < other:
                #    other, mrb2 = mrb2, other
                #
                mrb2 = mrb2 - 1.0
                temp = temp / ((mrb2 + mxib1) * mnib1)
                ccbnb = ccbnb + temp

            if mrb2 >= 1.0:
                return ccbnb
            temp = temp * mrb2 / (mnib1 + mrb2)
        else:
            mxib1 = beta_shape1
            mnib1 = i
            if temp == 0.0:
                temp = pmf_BetaNegativeBinomial(mnib1, mrb2, mxib1, other)
            else:
                temp = temp * mrb2 * other / (i * (mrb2 + other + mnib1 + mxib1 + -1))

            ccbnb = ccbnb + temp

        max_iterations = 60.0
        while True:
            temp = (
                temp * (mnib1 + mrb2) * (mnib1 + other) / (mnib1 + mxib1 + mrb2 + other)
            )
            mnib1 = mnib1 + 1.0
            if mxib1 < mnib1:
                mxib1, mnib1 = mnib1, mxib1

            temp = temp / mnib1
            ccbnb = ccbnb + temp
            if (temp <= 1e-16 * ccbnb) or (mnib1 + mxib1 > max_iterations):
                break

        temp = temp * (mnib1 + mrb2) * (mnib1 + other) / ((mnib1 + 1.0) * mxib1)
        ccbnb = ccbnb + temp
        mnib1 = mnib1 + 1.0
        mrb2 = mrb2 - 1.0
        other = other - 1.0
        while True:
            mnib1 = mnib1 + 1.0
            mxib1 = mxib1 + 1.0
            temp = temp * (mrb2 * other) / (mnib1 * mxib1)
            mrb2 = mrb2 - 1.0
            other = other - 1.0
            ccbnb = ccbnb + temp
            if abs(temp) <= 1e-16 * ccbnb:
                break

    return ccbnb


def pmf_BetaBinomial(i, sample_size, beta_shape1, beta_shape2):
    i = AlterForIntegralChecks_Others(i)
    sample_size = AlterForIntegralChecks_Others(sample_size)
    # These probably should be checked once, when creating a BetaBinomial distribution object
    if (beta_shape1 <= 0.0) or (beta_shape2 <= 0.0) or (sample_size < 0.0):
        raise ValueError
    elif (i < 0) or (i > sample_size):
        return 0.0
    else:
        return (
            (beta_shape1 / (i + beta_shape1))
            * (beta_shape2 / (beta_shape1 + beta_shape2))
            * (
                (sample_size + beta_shape1 + beta_shape2)
                / (sample_size - i + beta_shape2)
            )
            * hypergeometricTerm(i, sample_size - i, beta_shape1, beta_shape2)
        )


def cdf_BetaBinomial(i, sample_size, beta_shape1, beta_shape2):
    i = math.floor(i)
    sample_size = AlterForIntegralChecks_Others(sample_size)
    if (beta_shape1 <= 0.0) or (beta_shape2 <= 0.0) or (sample_size < 0.0):
        raise ValueError
    elif i < 0.0:
        return 0.0
    else:
        i = i + 1.0
        return sf_BetaNegativeBinomial(sample_size - i, i, beta_shape1, beta_shape2)


def sf_BetaBinomial(i, sample_size, beta_shape1, beta_shape2):
    i = math.floor(i)
    sample_size = AlterForIntegralChecks_Others(sample_size)
    if (beta_shape1 <= 0.0) or (beta_shape2 <= 0.0) or (sample_size < 0.0):
        raise ValueError
    elif i < 0.0:
        return 1.0
    elif i >= sample_size:
        return 0.0
    else:
        return sf_BetaNegativeBinomial(i, sample_size - i, beta_shape2, beta_shape1)


# Hypergeometric & Negative Hpergeometric pmfs, cdfs & sfs for completeness
def pmf_hypergeometric(type1s, sample_size, tot_type1, pop_size):
    type1s = AlterForIntegralChecks_Others(type1s)
    sample_size = AlterForIntegralChecks_Others(sample_size)
    tot_type1 = AlterForIntegralChecks_Others(tot_type1)
    pop_size = AlterForIntegralChecks_Others(pop_size)
    if (
        (sample_size < 0.0)
        or (tot_type1 < 0.0)
        or (sample_size > pop_size)
        or (tot_type1 > pop_size)
    ):
        raise ValueError
    else:
        return hypergeometricTerm(
            type1s,
            sample_size - type1s,
            tot_type1 - type1s,
            pop_size - tot_type1 - sample_size + type1s,
        )


def cdf_hypergeometric(type1s, sample_size, tot_type1, pop_size):
    type1s = math.floor(type1s)
    sample_size = AlterForIntegralChecks_Others(sample_size)
    tot_type1 = AlterForIntegralChecks_Others(tot_type1)
    pop_size = AlterForIntegralChecks_Others(pop_size)
    if (
        (sample_size < 0.0)
        or (tot_type1 < 0.0)
        or (sample_size > pop_size)
        or (tot_type1 > pop_size)
    ):
        raise ValueError
    else:
        return hypergeometric(
            type1s,
            sample_size - type1s,
            tot_type1 - type1s,
            pop_size - tot_type1 - sample_size + type1s,
            False,
        )


def sf_hypergeometric(type1s, sample_size, tot_type1, pop_size):
    type1s = math.floor(type1s)
    sample_size = AlterForIntegralChecks_Others(sample_size)
    tot_type1 = AlterForIntegralChecks_Others(tot_type1)
    pop_size = AlterForIntegralChecks_Others(pop_size)
    if (
        (sample_size < 0.0)
        or (tot_type1 < 0.0)
        or (sample_size > pop_size)
        or (tot_type1 > pop_size)
    ):
        raise ValueError
    else:
        return hypergeometric(
            type1s,
            sample_size - type1s,
            tot_type1 - type1s,
            pop_size - tot_type1 - sample_size + type1s,
            True,
        )


def pmf_neghypergeometric(type2s, type1s_reqd, tot_type1, pop_size):
    type2s = AlterForIntegralChecks_Others(type2s)
    type1s_reqd = AlterForIntegralChecks_Others(type1s_reqd)
    tot_type1 = AlterForIntegralChecks_Others(tot_type1)
    pop_size = AlterForIntegralChecks_Others(pop_size)
    if (type1s_reqd <= 0.0) or (tot_type1 < type1s_reqd) or (tot_type1 > pop_size):
        raise ValueError
    elif (type2s < 0.0) or (tot_type1 + type2s > pop_size):
        if type2s == 0.0:
            return 1.0
        else:
            return 0.0
    else:
        return (
            hypergeometricTerm(
                type1s_reqd - 1.0,
                type2s,
                tot_type1 - type1s_reqd + 1.0,
                pop_size - tot_type1 - type2s,
            )
            * (tot_type1 - type1s_reqd + 1.0)
            / (pop_size - type1s_reqd - type2s + 1.0)
        )


def cdf_neghypergeometric(type2s, type1s_reqd, tot_type1, pop_size):
    type2s = math.floor(type2s)
    type1s_reqd = AlterForIntegralChecks_Others(type1s_reqd)
    tot_type1 = AlterForIntegralChecks_Others(tot_type1)
    pop_size = AlterForIntegralChecks_Others(pop_size)
    if (type1s_reqd <= 0.0) or (tot_type1 < type1s_reqd) or (tot_type1 > pop_size):
        raise ValueError
    elif tot_type1 + type2s > pop_size:
        return 1.0
    else:
        return hypergeometric(
            type2s,
            type1s_reqd,
            pop_size - tot_type1 - type2s,
            tot_type1 - type1s_reqd,
            False,
        )


def sf_neghypergeometric(type2s, type1s_reqd, tot_type1, pop_size):
    type2s = math.floor(type2s)
    type1s_reqd = AlterForIntegralChecks_Others(type1s_reqd)
    tot_type1 = AlterForIntegralChecks_Others(tot_type1)
    pop_size = AlterForIntegralChecks_Others(pop_size)
    if (type1s_reqd <= 0.0) or (tot_type1 < type1s_reqd) or (tot_type1 > pop_size):
        raise ValueError
    elif tot_type1 + type2s > pop_size:
        return 0.0
    else:
        return hypergeometric(
            type2s,
            type1s_reqd,
            pop_size - tot_type1 - type2s,
            tot_type1 - type1s_reqd,
            True,
        )
