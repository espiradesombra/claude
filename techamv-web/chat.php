<?php
// Permitimos que la web lea la respuesta como JSON
header('Content-Type: application/json');

// Recogemos el texto que nos envía el frontend (tu web ZypyZape)
$datos = json_decode(file_get_contents('php://input'), true);
$textoUsuario = $datos['text'] ?? '';

if (empty($textoUsuario)) {
    echo json_encode(['reply' => 'No he recibido ningún mensaje.']);
    exit;
}

// === PON AQUÍ TU API KEY REAL DE GOOGLE AI STUDIO ===
$apiKey = 'AIzaSyBx4qHQEQlXkSZN3SbETgzG7OIoIToNHnM';
//'TU_API_KEY_REAL_AQUI';
// ====================================================

// La URL de Gemini (usando el modelo flash más rápido)
//$url = 'https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash-preview-09-2025:generateContent?key=' . $apiKey;
$url = 'https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key=' . $apiKey;
//$url = 'https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key=' . $apiKey;
// El comportamiento estricto del Asesor de ZypyZape
$systemPrompt = "Eres el portavoz y asesor oficial de ZypyZape, una revolucionaria tecnología de energía eólica basada en la 'Inercia de Fase'. 
      
TU OBJETIVO Y MENSAJES CLAVE:
- Explicar el salto de paradigma: ZypyZape ya no es un 'generador pasivo' que depende del capricho del viento, sino un 'regulador activo de red' que aporta estabilidad constante gracias a su inercia.
- Destacar el lucrativo 'Mercado de Regulación': Las energéticas ahora ganan dinero no solo vendiendo energía, sino vendiendo estabilidad y regulación de frecuencia a la red.
- La Jubilación Nuclear: Explicar que la energía nuclear se mantiene hoy porque aporta inercia base pesada para regular la red. ZypyZape aporta esta inercia de forma limpia y mecánica, haciendo que la nuclear quede obsoleta e innecesaria.
- Hablar del ahorro de gasto mundial (baterías, curtailment, condensadores síncronos, expansión de red).

REGLAS ESTRICTAS:
- NUNCA des detalles técnicos de ingeniería o mecánicos.
- NUNCA uses fórmulas matemáticas.
- Usa un tono visionario, corporativo, optimista y centrado en la ecología y el negocio (tonos claros, lenguaje ejecutivo verde).
- Sé conciso y claro en tus respuestas.";

// Preparamos el paquete de datos tal y como lo pide Google
$payload = [
    "contents" => [
        ["parts" => [["text" => $textoUsuario]]]
    ],
    "systemInstruction" => [
        "parts" => [["text" => $systemPrompt]]
    ]
];

// Iniciamos cURL para hacer la petición HTTP desde tu servidor en Nominalia
$ch = curl_init($url);
curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
curl_setopt($ch, CURLOPT_POST, true);
curl_setopt($ch, CURLOPT_HTTPHEADER, ['Content-Type: application/json']);
curl_setopt($ch, CURLOPT_POSTFIELDS, json_encode($payload));

// Ejecutamos y cerramos
$respuestaGoogle = curl_exec($ch);
curl_close($ch);

// Extraemos el texto de la respuesta de la IA
//$datosRespuesta = json_decode($respuestaGoogle, true);
//$respuestaBot = $datosRespuesta['candidates'][0]['content']['parts'][0]['text'] ?? "Lo siento, en este momento el servidor está //procesando mucha inercia de datos. Inténtalo en unos segundos.";
$datosRespuesta = json_decode($respuestaGoogle, true);

// Si Google nos devuelve un error, que nos lo diga por el chat
if (isset($datosRespuesta['error'])) {
    $respuestaBot = "ERROR DE GOOGLE: " . $datosRespuesta['error']['message'];
} else {
    $respuestaBot = $datosRespuesta['candidates'][0]['content']['parts'][0]['text'] ?? "Error desconocido. Respuesta de Google: " . $respuestaGoogle;
}
// Enviamos la respuesta de vuelta a la web web
echo json_encode(['reply' => $respuestaBot]);
?>