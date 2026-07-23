<?php
/**
 * FACTORITZACIÓ VMA - CORE SERVER-SIDE
 * Implementació del Mètode Diofàntic Cinemàtic (MDC) Autèntic
 * Basat en l'equació fonamental: N mod 2D = D
 */
header('Content-Type: application/json');
header("Access-Control-Allow-Origin: *");
header("Access-Control-Allow-Methods: POST");

// Rebre N (el semiprimer a factoritzar)
$data = json_decode(file_get_contents('php://input'), true);
$n_str = isset($data['n']) ? trim($data['n']) : '';

if (empty($n_str) || !ctype_digit($n_str)) {
    echo json_encode(['status' => 'error', 'message' => 'Cifra no vàlida']);
    exit;
}

// Forcem sencer de 64 bits estricte per màxima precisió (fins a 19 dígits)
$n = (int)$n_str; 

if ((string)$n !== $n_str) {
    echo json_encode(['status' => 'error', 'message' => 'La cifra supera el límit natiu de 64 bits del servidor.']);
    exit;
}

$start_time = microtime(true);

function factorize_mdc($num) {
    if ($num < 2) return null;
    // Purguem factors parells ràpidament
    if ($num % 2 === 0) return [2, intdiv($num, 2)];
    
    // DEFINICIÓ DELS LÍMITS MDC (Distància màxima = arrel de N)
    $root = (int)sqrt($num);
    $m_inicial = 1; 
    $m_final = intdiv($root - 3, 2); // m_max segons el document
    
    // BUCLE MDC AUTÈNTIC
    // Avalua l'espai L1 = {D = 2m+3 : m ∈ ℕ}
    for ($m = $m_inicial; $m <= $m_final; $m++) {
        $D = 2 * $m + 3; // Divisor candidat (imparell ≥ 5)
        
        // IDENTITAT DEL DOBLE MÒDUL VMA: N mod (2D) == D
        if ($num % (2 * $D) === $D) {
            return [$D, intdiv($num, $D)];
        }
        
        // Timeout preventiu del servidor (4.8 segons)
        if ((microtime(true) - $_SERVER['REQUEST_TIME_FLOAT']) > 4.8) {
            return "timeout";
        }
    }
    
    return null;
}

$result = factorize_mdc($n);

if ($result === "timeout") {
    echo json_encode(['status' => 'limit', 'message' => 'L\'espai de cerca requereix més temps de CPU (Timeout).']);
} elseif ($result) {
    $p = $result[0];
    $q = $result[1];
    
    $p_min = min($p, $q);
    $distancia = abs($p - $q);
    
    // Càlcul del límit de seguretat teòric (Salt Predictiu VMA)
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
        'method' => 'MDC_Core' // Afegit per demostrar a la terminal quin motor s'usa
    ]);
} else {
    echo json_encode(['status' => 'prime', 'message' => 'Cifra identificada com a Primer Atòmic.']);
}
?>