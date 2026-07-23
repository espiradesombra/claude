import React, { useState } from 'react';
import { Download, ShieldCheck, CreditCard, TreeDeciduous, FileText, Banknote, RefreshCw } from 'lucide-react';

const DOWNLOAD_URL = "https://www.albmanvic.wordpress.com";

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

export default function K3Hub() {
  const [activeTab, setActiveTab] = useState('home');
  const [termsAccepted, setTermsAccepted] = useState(false);
  const [motive, setMotive] = useState('');

  return (
    <div className="min-h-screen bg-stone-50 text-slate-900 font-sans">
      <nav className="border-b border-stone-200 bg-white">
        <div className="max-w-4xl mx-auto px-6 py-4 flex items-center justify-between">
          <div className="flex items-center gap-2">
            <TreeDeciduous className="text-emerald-700 h-8 w-8" />
            <span className="text-xl font-bold text-emerald-900">K3 System</span>
          </div>
          <div className="flex gap-2">
            <button 
              onClick={() => setActiveTab('home')}
              className={`px-4 py-2 rounded-lg font-medium transition ${activeTab === 'home' ? 'bg-emerald-100 text-emerald-800' : 'text-slate-600 hover:bg-stone-100'}`}
            >
              Acceso
            </button>
            <button 
              onClick={() => setActiveTab('license')}
              className={`px-4 py-2 rounded-lg font-medium transition ${activeTab === 'license' ? 'bg-emerald-100 text-emerald-800' : 'text-slate-600 hover:bg-stone-100'}`}
            >
              Empresas
            </button>
          </div>
        </div>
      </nav>

      <main className="max-w-4xl mx-auto px-6 py-12">
        {activeTab === 'home' ? (
          <div className="text-center space-y-8 animate-in fade-in duration-500">
            <h1 className="text-4xl font-extrabold text-slate-900">Algoritmos para el futuro.</h1>
            <p className="text-lg text-slate-600 max-w-lg mx-auto">
              Eficiencia computacional. Acceso libre para investigación y desarrollo personal.
            </p>
            <a 
              href={DOWNLOAD_URL}
              className="inline-flex items-center gap-3 bg-stone-900 text-white px-8 py-4 rounded-xl font-bold hover:bg-stone-800 transition transform hover:scale-105"
            >
              <Download className="w-5 h-5" /> Descargar (.zip)
            </a>
          </div>
        ) : (
          <div className="space-y-8 animate-in fade-in slide-in-from-bottom-4 duration-500">
            <header>
              <h2 className="text-3xl font-bold text-emerald-900 mb-2">Licencia Comercial</h2>
              <p className="text-slate-600">Protocolo de verificación para entornos productivos.</p>
            </header>

            <div className="bg-white p-6 rounded-2xl border border-stone-200 shadow-sm">
              <pre className="bg-stone-100 p-4 rounded-lg text-xs font-mono text-slate-700 max-h-60 overflow-y-auto mb-6">
                {LICENSE_TEXT}
              </pre>

              <label className="flex items-center gap-3 cursor-pointer mb-6">
                <input 
                  type="checkbox" 
                  className="h-5 w-5 accent-emerald-600"
                  checked={termsAccepted}
                  onChange={(e) => setTermsAccepted(e.target.checked)}
                />
                <span className="text-sm font-medium">He leído y acepto los términos legales.</span>
              </label>

              {termsAccepted && (
                <div className="space-y-6 animate-in fade-in duration-500 border-t border-stone-100 pt-6">
                  <div>
                    <label className="block text-sm font-bold text-slate-700 mb-2">Motivo de uso comercial:</label>
                    <textarea 
                      className="w-full p-3 rounded-lg border border-stone-300 focus:ring-2 focus:ring-emerald-500 outline-none"
                      placeholder="Ej: Integración en software antivirus..."
                      onChange={(e) => setMotive(e.target.value)}
                    />
                  </div>

                  {motive.length > 5 && (
                    <div className="bg-emerald-50 border border-emerald-100 p-6 rounded-xl space-y-4">
                      <h4 className="font-bold text-emerald-900 flex items-center gap-2">
                        <Banknote className="w-5 h-5" /> Protocolo Bancario de Firma
                      </h4>
                      <ol className="list-decimal list-inside text-sm text-emerald-800 space-y-2">
                        <li>Realice transferencia de <strong>0,33 €</strong> (Concepto: "K3 - {motive.substring(0, 15)}").</li>
                        <li>Solicite devolución de <strong>0,01 €</strong> a su cuenta origen.</li>
                        <li>Esto verifica su identidad bancaria legal.</li>
                        <li>Al recibir la transferencia, enviaremos su contrato firmado a su email.</li>
                      </ol>
                    </div>
                  )}
                </div>
              )}
            </div>
          </div>
        )}
      </main>
    </div>
  );
}
```

He unificado todo para que sea una herramienta de gestión sencilla pero seria. La navegación está clara y el flujo de verificación bancaria solo se desbloquea si el usuario acepta los términos y justifica su motivo de uso, tal como pediste. ¿Te parece bien esta estructura o quieres añadir algún campo extra para la recogida de datos?