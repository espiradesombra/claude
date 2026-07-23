<?php
header('Content-Type: application/json');
header("Access-Control-Allow-Origin: *");
header("Access-Control-Allow-Methods: POST");

$data = json_decode(file_get_contents('php://input'), true);
$n_str = isset($data['n']) ? trim($data['n']) : '';

if (empty($n_str) || !ctype_digit($n_str)) {
    echo json_encode(['status' => 'error', 'message' => 'Número inválido']);
    exit;
}

$n = (int)$n_str;
if ((string)$n !== $n_str) {
    echo json_encode(['status' => 'error', 'message' => 'Overflow 64-bit']);
    exit;
}

$start = microtime(true);

/* =========================
   CONFIG
========================= */
$PRIMOS = [2,3,5,7,11,13,17,19,23];
$PRIMORIAL = array_product($PRIMOS);
$LIM_CRIBA = 2000000;
$MAX_BISEC = 8000;
$TALL_ASIM = 4;

/* =========================
   UTILIDADES
========================= */

function pasa_rueda($m, $primorial) {
    if ($m < 1) return false;
    return gcd(2*$m+3, $primorial) === 1;
}

function gcd($a, $b) {
    while ($b != 0) {
        $t = $b;
        $b = $a % $b;
        $a = $t;
    }
    return $a;
}

function proper_valid($m, $pas, $primorial) {
    for ($i=0; $i<200; $i++) {
        if ($m < 1) return null;
        if (pasa_rueda($m, $primorial)) return $m;
        $m += $pas;
    }
    return null;
}

/* =========================
   ΔΦ SIGNO (núcleo)
========================= */

function sgn_desfase($m, $N) {
    $B = 2*$m + 3;
    if ($B <= 1 || $B >= $N) return 0;

    $D_B = 2*$B;
    $R_B = $N % $D_B;

    $A = intdiv($N, $B);
    if ($A < 3) return 0;

    $m_e = intdiv($A - 3, 2);
    if ($m_e < 1) return 0;

    $D_A = 2*(2*$m_e + 3);
    $R_A = $N % $D_A;

    $d = $R_B * $D_A - $R_A * $D_B;

    if ($d > 0) return 1;
    if ($d < 0) return -1;
    return 0;
}

function check_factor($m, $N) {
    $D = 2*$m + 3;
    return ($D > 1 && $D < $N && $N % $D === 0);
}

/* =========================
   BISECCIÓN ASIMÉTRICA
========================= */

function biseccion_fase($N, $m_ini, $m_fi, $primorial) {
    global $MAX_BISEC, $TALL_ASIM;

    for ($it=0; $it < $MAX_BISEC; $it++) {

        $rang = $m_fi - $m_ini;
        if ($rang < 50) {
            // fallback directo corto
            for ($m=$m_ini; $m<=$m_fi; $m++) {
                if (check_factor($m, $N)) return 2*$m+3;
            }
            return null;
        }

        $m_25 = proper_valid($m_ini + intdiv($rang, $TALL_ASIM), 1, $primorial);
        $m_e  = proper_valid($m_ini, 1, $primorial);
        $m_d  = proper_valid($m_fi, -1, $primorial);

        if (!$m_25 || !$m_e || !$m_d) return null;

        // chequeo directo
        foreach ([$m_e, $m_25, $m_d] as $mp) {
            if (check_factor($mp, $N)) return 2*$mp+3;
        }

        $s_e  = sgn_desfase($m_e,  $N);
        $s_25 = sgn_desfase($m_25, $N);
        $s_d  = sgn_desfase($m_d,  $N);

        // decisión
        if ($s_e != 0 && $s_25 != 0 && $s_e != $s_25) {
            $m_fi = $m_25; // zona 25%
        }
        elseif ($s_25 != 0 && $s_d != 0 && $s_25 != $s_d) {
            $m_ini = $m_25; // zona 75%
        }
        else {
            // sin señal → probablemente cerca de sqrt o plano
            return null;
        }

        // timeout duro
        if ((microtime(true) - $_SERVER['REQUEST_TIME_FLOAT']) > 4.8) {
            return "timeout";
        }
    }

    return null;
}

/* =========================
   ORQUESTADOR
========================= */

function factor_mdc_fase($N) {
    global $PRIMOS, $PRIMORIAL, $LIM_CRIBA;

    // F0 — primos pequeños
    foreach ($PRIMOS as $p) {
        if ($N % $p === 0 && $p < $N) return $p;
    }

    $r = (int)sqrt($N);
    if ($r*$r === $N) return $r;

    $m_max = intdiv($r - 3, 2);
    $lim = min($m_max, $LIM_CRIBA);

    // criba
    for ($m=1; $m <= $lim; $m++) {
        if (pasa_rueda($m, $PRIMORIAL)) {
            $D = 2*$m+3;
            if ($N % $D === 0) return $D;
        }
    }

    if ($lim >= $m_max) return null;

    // FASE PURA (sin K-sweep)
    return biseccion_fase($N, $lim+1, $m_max, $PRIMORIAL);
}

/* =========================
   EJECUCIÓN
========================= */

$res = factor_mdc_fase($n);

if ($res === "timeout") {
    echo json_encode(['status'=>'limit']);
}
elseif ($res) {
    echo json_encode([
        'status'=>'success',
        'p'=>(string)$res,
        'q'=>(string)intdiv($n,$res),
        'method'=>'MDC_FASE_PURA'
    ]);
}
else {
    echo json_encode([
        'status'=>'unknown',
        'message'=>'No detectado por fase (posible zona √N o señal plana)'
    ]);
}
?>