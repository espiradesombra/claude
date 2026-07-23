<?php
header('Content-Type: application/json');

$data = json_decode(file_get_contents('php://input'), true);
$n_str = isset($data['n']) ? trim($data['n']) : '';

if (!ctype_digit($n_str)) {
    echo json_encode(['status'=>'error']); exit;
}

$n = (int)$n_str;
if ((string)$n !== $n_str) {
    echo json_encode(['status'=>'overflow']); exit;
}

/* ================= CONFIG ================= */

$PRIMOS = [2,3,5,7,11,13,17,19,23];
$PRIMORIAL = array_product($PRIMOS);
$LIM_CRIBA = 2000000;
$MAX_IT = 4000;

/* ================= UTIL ================= */

function gcd($a,$b){while($b){$t=$b;$b=$a%$b;$a=$t;}return $a;}

function pasa_rueda($m,$primorial){
    return ($m>=1 && gcd(2*$m+3,$primorial)===1);
}

function check_factor($m,$N){
    $D = 2*$m+3;
    return ($D>1 && $D<$N && $N % $D === 0);
}

/* ================= |ΔΦ| ================= */

function mag_desfase($m,$N){
    $B = 2*$m + 3;
    if ($B <= 1 || $B >= $N) return PHP_INT_MAX;

    $D_B = 2*$B;
    $R_B = $N % $D_B;

    $A = intdiv($N,$B);
    if ($A < 3) return PHP_INT_MAX;

    $m_e = intdiv($A - 3, 2);
    if ($m_e < 1) return PHP_INT_MAX;

    $D_A = 2*(2*$m_e + 3);
    $R_A = $N % $D_A;

    return abs($R_B * $D_A - $R_A * $D_B);
}

/* ================= BÚSQUEDA POR MÍNIMO ================= */

function buscar_minimo($N,$m_ini,$m_fi,$primorial){

    global $MAX_IT;

    for ($it=0;$it<$MAX_IT;$it++){

        $rang = $m_fi - $m_ini;
        if ($rang < 50){
            for($m=$m_ini;$m<=$m_fi;$m++){
                if(check_factor($m,$N)) return 2*$m+3;
            }
            return null;
        }

        // muestreo 3 puntos (puedes subir a 5 si quieres más precisión)
        $m1 = $m_ini;
        $m2 = intdiv($m_ini+$m_fi,2);
        $m3 = $m_fi;

        $candidatos = [];

        foreach([$m1,$m2,$m3] as $m){
            if(!pasa_rueda($m,$primorial)) continue;

            if(check_factor($m,$N)) return 2*$m+3;

            $val = mag_desfase($m,$N);
            $candidatos[] = [$m,$val];
        }

        if(empty($candidatos)) return null;

        // elegir mínimo
        usort($candidatos,function($a,$b){
            return $a[1] <=> $b[1];
        });

        $best_m = $candidatos[0][0];

        // contraer rango alrededor del mínimo
        $delta = intdiv($rang,4); // contracción fuerte

        $m_ini = max(1, $best_m - $delta);
        $m_fi  = $best_m + $delta;

        // timeout
        if ((microtime(true) - $_SERVER['REQUEST_TIME_FLOAT']) > 4.8){
            return "timeout";
        }
    }

    return null;
}

/* ================= CORE ================= */

function factor_mdc_min($N){

    global $PRIMOS,$PRIMORIAL,$LIM_CRIBA;

    foreach($PRIMOS as $p){
        if($N % $p === 0 && $p < $N) return $p;
    }

    $r = (int)sqrt($N);
    if($r*$r === $N) return $r;

    $m_max = intdiv($r-3,2);
    $lim = min($m_max,$LIM_CRIBA);

    // criba
    for($m=1;$m<=$lim;$m++){
        if(pasa_rueda($m,$PRIMORIAL)){
            $D = 2*$m+3;
            if($N % $D === 0) return $D;
        }
    }

    if($lim >= $m_max) return null;

    return buscar_minimo($N,$lim+1,$m_max,$PRIMORIAL);
}

/* ================= RUN ================= */

$res = factor_mdc_min($n);

if ($res === "timeout"){
    echo json_encode(['status'=>'limit']);
}
elseif ($res){
    echo json_encode([
        'status'=>'success',
        'p'=>(string)$res,
        'q'=>(string)intdiv($n,$res),
        'method'=>'MDC_MIN_DELTA'
    ]);
}
else{
    echo json_encode([
        'status'=>'fail',
        'message'=>'mínimo no suficientemente definido'
    ]);
}
?>