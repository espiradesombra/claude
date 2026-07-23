<?php
// Permitimos que la web lea la respuesta como JSON
header('Content-Type: application/json');
header("Access-Control-Allow-Origin: *");
header("Access-Control-Allow-Methods: POST");

// Recogemos el texto que nos envía el formulario de React
$datos = json_decode(file_get_contents('php://input'), true);

$nombre = $datos['nombre'] ?? '';
$email = $datos['email'] ?? '';
$empresa = $datos['empresa'] ?? '';
$pais = $datos['pais'] ?? '';
$intencion = $datos['intencion'] ?? '';
$mensaje = $datos['mensaje'] ?? '';

// Validación básica (comprobar que los campos obligatorios no están vacíos)
if (empty($nombre) || empty($email) || empty($mensaje) || empty($intencion)) {
    http_response_code(400);
    echo json_encode(['error' => 'Faltan campos obligatorios.']);
    exit;
}

// === CONFIGURACIÓN DE LOS CORREOS ===
$destinatarios = "victor@techamv.com, david@techamv.com";
$asunto = "🚀 Nuevo Contacto ZypyZape (Web) - Interés: " . strtoupper($intencion);

// Construimos el cuerpo del correo de forma elegante
$cuerpo = "Has recibido una nueva solicitud de contacto desde la web de ZypyZape.\n\n";
$cuerpo .= "======================================\n";
$cuerpo .= "👤 Nombre / Entidad: " . $nombre . "\n";
$cuerpo .= "📧 Email de contacto: " . $email . "\n";
$cuerpo .= "🏢 Empresa / Sector: " . ($empresa ? $empresa : 'No especificado') . "\n";
$cuerpo .= "🌍 País: " . ($pais ? $pais : 'No especificado') . "\n";
$cuerpo .= "🎯 Intención: " . $intencion . "\n";
$cuerpo .= "======================================\n\n";
$cuerpo .= "📝 DETALLES DEL MENSAJE:\n";
$cuerpo .= $mensaje . "\n\n";
$cuerpo .= "======================================\n";
$cuerpo .= "Mensaje automático enviado desde techamv.com";

// Cabeceras del correo (Importante para que no vaya a Spam)
// Se enviará a través del servidor local, y si el usuario le da a "Responder", os responderá al email del cliente
$headers = "From: web@techamv.com\r\n"; 
$headers .= "Reply-To: " . $email . "\r\n";
$headers .= "X-Mailer: PHP/" . phpversion();

// Intentamos enviar el correo mediante la función nativa de PHP
$enviado = mail($destinatarios, $asunto, $cuerpo, $headers);

// Respondemos a la web si ha sido un éxito o un error
if ($enviado) {
    echo json_encode(['success' => true, 'message' => 'Correo enviado correctamente.']);
} else {
    http_response_code(500);
    echo json_encode(['error' => 'Error al enviar el correo. Revisa la función mail() del hosting.']);
}
?>