import React, { useState } from 'react';
import { Leaf, ShieldCheck, Download, CreditCard, Code, TreeDeciduous, CheckSquare, FileText, ArrowRight, Lock } from 'lucide-react';

const LICENSE_TEXT = `LICENCIA DE USO DE SOFTWARE K3

Este es un acuerdo legal entre usted (el "Licenciatario") y Víctor (el "Licenciante") para el uso de la tecnología, algoritmos y binarios asociados al proyecto "K3".

### 1. ACEPTACIÓN DE TÉRMINOS
Al realizar el pago de la tarifa establecida (133 €) o al acceder al enlace de uso gratuito, usted acepta expresamente los términos de esta Licencia. El pago de 133 € constituye su firma digital y aceptación vinculante de este contrato.

### 2. MODALIDADES DE LICENCIA
A. Licencia de Uso Personal (Gratuita):
Se concede una licencia limitada, no exclusiva, gratuita y para uso exclusivamente personal, académico o de investigación no comercial.
* Usted puede hashear, probar y usar las herramientas para proyectos personales.
* Queda prohibido su uso en entornos productivos, comerciales, o integrados en software de terceros sin la licencia de pago.

B. Licencia de Uso Comercial (133 € - De por vida):
Para cualquier uso con fines de lucro, comercial, o integración en productos de terceros (ejemplos: antivirus, motores de seguridad, sistemas de cifrado, análisis de datos comerciales), es obligatorio adquirir la Licencia Comercial.
* Tarifa: 133 € (Pago único).
* Vigencia: Perpetua (de por vida).
* Derechos: Autoriza el uso del código y binarios en entornos comerciales y empresariales.
* Incumplimiento: El uso comercial de la tecnología sin la Licencia Comercial activa constituye una violación de los derechos de autor y propiedad intelectual del Licenciante.

### 3. PROPIEDAD INTELECTUAL Y GARANTÍAS
El software se proporciona "tal cual" (AS-IS). El Licenciante no ofrece garantías de funcionamiento para fines específicos y no se hace responsable de daños derivados del uso del software. La propiedad intelectual del algoritmo y la teoría subyacente pertenecen exclusivamente al Licenciante.

### 4. JURISDICCIÓN
Este acuerdo se rige por las leyes vigentes en España. Cualquier disputa será resuelta en los tribunales correspondientes al domicilio del Licenciante.
Al proceder con el pago, usted declara ser mayor de edad y haber leído y comprendido los términos anteriores.`;

export default function App() {
  const [termsAccepted, setTermsAccepted] = useState(false);
  const [showVerification, setShowVerification] = useState(false);
  const [motive, setMotive] = useState('');

  return (
    <div className="min-h-screen bg-stone-50 text-slate-900 font-sans selection:bg-emerald-200">
      
      <nav className="p-6 flex justify-between items-center max-w-6xl mx-auto">
        <div className="flex items-center gap-2">
          <TreeDeciduous className="text-emerald-700 h-8 w-8" />
          <span className="text-2xl font-bold tracking-tight text-emerald-900">K3 Software</span>
        </div>
      </nav>

      <main className="max-w-4xl mx-auto px-6 py-12">
        <header className="text-center mb-16">
          <h1 className="text-4xl md:text-5xl font-extrabold text-slate-900 mb-6">
            Algoritmos <span className="text-emerald-700">inteligentes</span>.
          </h1>
          <p className="text-lg text-slate-600 max-w-2xl mx-auto">
            Eficiencia computacional para un mundo sostenible. Elige tu modo de acceso.
          </p>
        </header>

        {}
        <div className="grid md:grid-cols-2 gap-8">
          
          {/* Free Tier */}
          <div className="bg-white p-8 rounded-2xl border border-stone-200 shadow-sm flex flex-col justify-between">
            <div>
              <h3 className="text-xl font-bold mb-4">Uso Personal</h3>
              <p className="text-slate-600 mb-6 text-sm">Ideal para investigadores y proyectos personales. Acceso sin registro.</p>
            </div>
            <button className="w-full flex items-center justify-center gap-2 bg-stone-900 hover:bg-stone-800 text-white py-3 rounded-lg font-bold transition">
              <Download className="w-4 h-4" /> Descargar (.zip)
            </button>
          </div>

          {/* Commercial Tier */}
          <div className="bg-emerald-900 p-8 rounded-2xl text-white shadow-xl flex flex-col justify-between relative overflow-hidden">
             <div className="absolute top-0 right-0 p-4 bg-emerald-800 rounded-bl-xl text-xs font-bold uppercase tracking-widest text-emerald-100">
               Certificado
             </div>
             <div>
              <h3 className="text-xl font-bold mb-4">Licencia Comercial</h3>
              <p className="text-emerald-100 mb-6 text-sm">Para entornos productivos y software de terceros. Requiere verificación bancaria.</p>
            </div>
            <button 
              onClick={() => setShowVerification(true)}
              className="w-full flex items-center justify-center gap-2 bg-emerald-500 hover:bg-emerald-400 text-emerald-950 py-3 rounded-lg font-bold transition"
            >
              <CreditCard className="w-4 h-4" /> Iniciar Verificación
            </button>
          </div>
        </div>

        {}
        {showVerification && (
          <div className="mt-12 bg-white p-8 rounded-2xl border border-emerald-200 shadow-2xl animate-in fade-in slide-in-from-bottom-4 duration-500">
            <h2 className="text-2xl font-bold mb-6 text-emerald-900 flex items-center gap-2">
              <ShieldCheck /> Protocolo de Verificación K3
            </h2>
            
            <div className="prose prose-stone max-w-none text-sm text-slate-600 mb-8">
              <pre className="bg-stone-100 p-4 rounded text-xs font-mono max-h-40 overflow-y-auto border border-stone-200">
                {LICENSE_TEXT}
              </pre>
            </div>

            <div className="space-y-6">
              <label className="flex items-start gap-3 cursor-pointer">
                <input 
                  type="checkbox" 
                  className="mt-1 accent-emerald-600 h-5 w-5"
                  checked={termsAccepted}
                  onChange={(e) => setTermsAccepted(e.target.checked)}
                />
                <span className="text-sm font-medium text-slate-800">
                  He leído, comprendido y acepto las condiciones legales descritas anteriormente.
                </span>
              </label>

              {termsAccepted && (
                <div className="space-y-4 animate-in fade-in duration-500">
                  <div className="bg-emerald-50 p-6 rounded-xl border border-emerald-100">
                    <h4 className="font-bold text-emerald-900 mb-2">1. Motivo de uso comercial:</h4>
                    <textarea 
                      className="w-full p-3 rounded-lg border border-emerald-200 focus:ring-2 focus:ring-emerald-500 outline-none"
                      placeholder="Ej: Integración en software antivirus para validación de ficheros..."
                      value={motive}
                      onChange={(e) => setMotive(e.target.value)}
                    />
                  </div>

                  <div className="bg-amber-50 p-6 rounded-xl border border-amber-100 text-amber-900">
                    <h4 className="font-bold mb-2">2. Instrucciones de Pago y Firma:</h4>
                    <ol className="list-decimal list-inside space-y-2 text-sm">
                      <li>Realice una transferencia de <strong>0,33 €</strong> a [IBAN-TU-REVOLUT].</li>
                      <li>Indique en concepto: "K3 - {motive.substring(0, 10) || 'Licencia'}"</li>
                      <li>Solicite una devolución de <strong>0,01 €</strong> (esto verifica su identidad bancaria).</li>
                      <li>Al recibir el pago, enviaremos a su email el contrato oficial para firma digital.</li>
                    </ol>
                  </div>
                </div>
              )}
            </div>
          </div>
        )}
      </main>

      <footer className="py-12 text-center text-stone-400 text-xs">
        &copy; {new Date().getFullYear()} Víctor. Todos los derechos reservados.
      </footer>
    </div>
  );
}