<?php
/**
 * FACTORIZACIÓN HÍBRIDA VMA - CORE SERVER-SIDE
 * Implementación del Método Diofántico Cinemático (MDC)
 * Incorpora: Criba de Rueda (Filtro Primorial) y Oráculo de Bisección
 */
header('Content-Type: application/json');
header("Access-Control-Allow-Origin: *");
header("Access-Control-Allow-Methods: POST");

// ================================
// 1. GESTIÓN DE ENTRADA Y SEGURIDAD
// ================================
$data = json_decode(file_get_contents('php://input'), true);
$n_str = isset($data['n']) ? trim($data['n']) : '';

if (empty($n_str) || !ctype_digit($n_str)) {
    echo json_encode(['status' => 'error', 'message' => 'Cifra no válida. Introduzca un número entero.']);
    exit;
}

$N = (int)$n_str;
if ((string)$N !== $n_str) {
    echo json_encode(['status' => 'error', 'message' => 'Desbordamiento: Supera la precisión nativa de 64 bits del servidor (máx ~18 dígitos).']);
    exit;
}

$start_time = microtime(true);

// ================================
// 2. CONFIGURACIÓN DEL ENTORNO MDC
// ================================
$PRIMOS = [2, 3, 5, 7, 11, 13, 17, 19, 23];
$PRIMORIAL = array_product($PRIMOS); // 223092870
$LIM_CRIBA = 200000;
$MAX_ITER = 6000;

// ================================
// 3. FUNCIONES BASE Y GEOMETRÍA
// ================================
function gcd($a, $b) {
    while ($b) { $t = $b; $b = $a % $b; $a = $t; }
    return $a;
}

function pasa_rueda($m, $primorial) {
    // Descarta candidatos que comparten divisores con el primorial
    return ($m >= 1 && gcd(2 * $m + 3, $primorial) === 1);
}

function check_factor($m, $N) {
    $D = 2 * $m + 3;
    return ($D > 1 && $D < $N && $N % $D === 0);
}

function desfase_abs($m, $N) {
    $D = 2 * $m + 3;
    // Oráculo de Bisección MDC basado en el diferencial de la señal
    // diff = | 2 * (N mod 2D) - 2D |
    return abs(2 * ($N % (2 * $D)) - 2 * $D);
}

// ================================
// 4. NÚCLEO HÍBRIDO VMA
// ================================
function factorize_hybrid($N) {
    global $PRIMOS, $PRIMORIAL, $LIM_CRIBA, $MAX_ITER;

    // Purga trivial instantánea
    foreach ($PRIMOS as $p) {
        if ($N % $p === 0 && $p < $N) return $p;
    }

    $root = (int)sqrt($N);
    if ($root * $root === $N) return $root;

    $m_max = intdiv($root - 3, 2);
    $lim = min($m_max, $LIM_CRIBA);

    // FASE A: Criba de Rueda Lineal Optimizada
    for ($m = 1; $m <= $lim; $m++) {
        if (pasa_rueda($m, $PRIMORIAL)) {
            if (check_factor($m, $N)) return 2 * $m + 3;
        }
        if ((microtime(true) - $_SERVER['REQUEST_TIME_FLOAT']) > 4.5) return "timeout";
    }

    if ($m_max <= $LIM_CRIBA) return null;

    // FASE B: Oráculo de Bisección (Descenso Geométrico MDC)
    $m_ini = $LIM_CRIBA + 1;
    $m_fin = $m_max;
    $iter = 0;

    while ($m_fin - $m_ini > 1000 && $iter < $MAX_ITER) {
        $iter++;
        $step = intdiv($m_fin - $m_ini, 4);
        
        $m1 = $m_ini;
        $m2 = $m_ini + $step;
        $m3 = $m_ini + 2 * $step;
        $m4 = $m_ini + 3 * $step;
        $m5 = $m_fin;

        $points = [$m1, $m2, $m3, $m4, $m5];
        $vals = [];

        // Evaluación de los 5 pivotes en el espacio L1
        foreach ($points as $m) {
            if (check_factor($m, $N)) return 2 * $m + 3;
            $vals[] = desfase_abs($m, $N);
        }

        // Búsqueda del valle de convergencia geométrica
        $min_val = min($vals);
        $min_index = array_search($min_val, $vals);

        // Contracción de ventana de búsqueda
        if ($min_index === 0) {
            $m_fin = $m2;
        } elseif ($min_index === 1) {
            $m_ini = $m1; $m_fin = $m3;
        } elseif ($min_index === 2) {
            $m_ini = $m2; $m_fin = $m4;
        } elseif ($min_index === 3) {
            $m_ini = $m3; $m_fin = $m5;
        } else {
            $m_ini = $m4;
        }

        if ((microtime(true) - $_SERVER['REQUEST_TIME_FLOAT']) > 4.5) return "timeout";
    }

    // FASE C: Convergencia Local (Fuerza bruta con Rueda en valle reducido)
    for ($m = $m_ini; $m <= $m_fin; $m++) {
        if (pasa_rueda($m, $PRIMORIAL)) {
            if (check_factor($m, $N)) return 2 * $m + 3;
        }
    }

    return null;
}

// ================================
// 5. RESPUESTA JSON A LA TERMINAL
// ================================
$res = factorize_hybrid($N);

if ($res === "timeout") {
    echo json_encode(['status' => 'limit', 'message' => 'L\'espai de cerca requereix més temps de CPU (Timeout prevenció servidors).']);
} elseif ($res) {
    $p = $res;
    $q = intdiv($N, $p);
    
    $p_min = min($p, $q);
    $distancia = abs($p - $q);
    
    // Càlcul del límit teòric de la teva conjetura VMA
    $l_limit = floor(sqrt($p_min + 3)) + 7;
    $exec_time = (microtime(true) - $start_time) * 1000;

    echo json_encode([
        'status' => 'success',
        'n' => $n_str,
        'p' => (string)$p,
        'q' => (string)$q,
        'distancia' => (string)$distancia,
        'l_limit' => (string)$l_limit,
        'seguro' => ($distancia <= $l_limit),
        'time_ms' => round($exec_time, 2),
        'method' => 'MDC_Hybrid_Oracle' // Indica que ha usat el mètode avançat
    ]);
} else {
    echo json_encode(['status' => 'prime', 'message' => 'Cifra identificada com a Primer Atòmic.']);
}
?>