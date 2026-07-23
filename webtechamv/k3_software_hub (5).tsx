import React, { useState, useEffect } from 'react';
import { Download, ShieldCheck, CreditCard, TreeDeciduous, FileText, Banknote, RefreshCw, Lock, CheckCircle } from 'lucide-react';

const DOWNLOAD_URL = "https://www.albmanvic.wordpress.com";

const CONTRACT_TEXT = `CONTRATO DE LICENCIA DE USO Y EXPLOTACIÓN VITALICIA

En [CIUDAD], a fecha actual.

DE UNA PARTE:
D. VÍCTOR [APELLIDOS], en adelante el LICENCIANTE, con domicilio a efectos de notificaciones en el teléfono +674 941 533.

DE OTRA PARTE:
[NOMBRE DE LA EMPRESA O RAZÓN SOCIAL], con NIF/CIF [NÚMERO], en adelante el LICENCIATARIO.

EXPONEN:
Que el Licenciante es titular de los derechos de explotación sobre el algoritmo "HashCode" aplicable al ámbito de la identificación.
Que el Licenciatario desea adquirir una licencia vitalicia.

CLÁUSULAS:
PRIMERA. OBJETO: Licencia no exclusiva, intransferible y vitalicia para el uso comercial de la tecnología "HashCode".
SEGUNDA. PRECIO: Pago único de 133,00 € (ciento treinta y tres euros) mediante transferencia/Bizum al +674 941 533.
TERCERA. PROPIEDAD INTELECTUAL: Propiedad total del Licenciante.
CUARTA. EXCLUSIÓN DE GARANTÍAS: Software "AS-IS".
QUINTA. JURISDICCIÓN: Tribunales de España.

FIRMAS:
EL LICENCIANTE: Amor
EL LICENCIATARIO: [Nombre y Cargo]`;

export default function K3Hub() {
  const [activeTab, setActiveTab] = useState('home');
  const [isVerifying, setIsVerifying] = useState(false);
  const [hashValue, setHashValue] = useState('');
  const [showContact, setShowContact] = useState(false);

  // Simulate hash generation for the "handshake" effect
  const handleVerify = () => {
    setIsVerifying(true);
    setShowContact(false);
    
    // Simulate technical processing
    setTimeout(() => {
      const mockHash = "0x" + Array.from({length: 64}, () => Math.floor(Math.random()*16).toString(16)).join('');
      setHashValue(mockHash);
      setIsVerifying(false);
      setShowContact(true);
    }, 2500);
  };

  return (
    <div className="min-h-screen bg-stone-50 text-slate-900 font-sans">
      <nav className="border-b border-stone-200 bg-white">
        <div className="max-w-4xl mx-auto px-6 py-4 flex items-center justify-between">
          <div className="flex items-center gap-2">
            <TreeDeciduous className="text-emerald-700 h-8 w-8" />
            <span className="text-xl font-bold text-emerald-900 tracking-tight">K3 System</span>
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
              Licencia Empresas
            </button>
          </div>
        </div>
      </nav>

      <main className="max-w-4xl mx-auto px-6 py-12">
        {activeTab === 'home' ? (
          <div className="text-center space-y-8 animate-in fade-in duration-500">
            <h1 className="text-5xl font-extrabold text-slate-900 tracking-tighter">Algoritmos para el futuro.</h1>
            <p className="text-xl text-slate-600 max-w-lg mx-auto">
              Eficiencia computacional. Acceso libre para investigación. Protección mediante HashCode.
            </p>
            <a 
              href={DOWNLOAD_URL}
              className="inline-flex items-center gap-3 bg-stone-900 text-white px-8 py-4 rounded-xl font-bold hover:bg-emerald-800 transition transform hover:scale-105 shadow-lg"
            >
              <Download className="w-5 h-5" /> Descargar (.zip)
            </a>
          </div>
        ) : (
          <div className="space-y-8 animate-in fade-in slide-in-from-bottom-4 duration-500">
            <header>
              <h2 className="text-3xl font-bold text-emerald-900 mb-2">Protocolo de Licencia Comercial</h2>
              <p className="text-slate-600">Verificación de identidad y validación de contrato.</p>
            </header>

            <div className="bg-white p-6 rounded-2xl border border-stone-200 shadow-sm space-y-6">
              <h3 className="font-bold text-lg text-slate-800">Contrato de Licencia</h3>
              <pre className="bg-stone-900 text-emerald-400 p-6 rounded-lg text-xs font-mono max-h-60 overflow-y-auto whitespace-pre-wrap border border-stone-800">
                {CONTRACT_TEXT}
              </pre>

              <div className="bg-emerald-50 p-6 rounded-xl border border-emerald-100 space-y-4">
                <h4 className="font-bold text-emerald-900 flex items-center gap-2">
                  <Lock className="w-5 h-5" /> Protocolo de Verificación de Identidad
                </h4>
                <p className="text-sm text-emerald-800">
                  Para proceder con la adquisición de la licencia, iniciamos un Handshake de seguridad bancaria.
                </p>

                {!showContact ? (
                  <button 
                    onClick={handleVerify}
                    disabled={isVerifying}
                    className="flex items-center gap-2 bg-emerald-600 text-white px-6 py-3 rounded-lg font-bold hover:bg-emerald-700 transition disabled:opacity-50"
                  >
                    {isVerifying ? <RefreshCw className="animate-spin" /> : <ShieldCheck />}
                    {isVerifying ? "Generando Handshake de Seguridad..." : "Iniciar Verificación (0,33€)"}
                  </button>
                ) : (
                  <div className="space-y-4 animate-in fade-in">
                    <div className="p-4 bg-emerald-900 rounded-lg font-mono text-xs text-emerald-300 break-all border border-emerald-700">
                      <span className="text-emerald-500 block mb-1">// Handshake Completo:</span>
                      {hashValue}
                    </div>
                    <div className="text-sm text-emerald-900 space-y-2">
                      <p><strong>Estado:</strong> Validado.</p>
                      <p>Realice la transferencia de <strong>0,33 €</strong> vía Bizum al <strong>+674 941 533</strong> para finalizar.</p>
                      <p>Tras la recepción, recibirá el contrato firmado por <strong>Amor</strong> en su correo.</p>
                    </div>
                  </div>
                )}
              </div>
            </div>
          </div>
        )}
      </main>
    </div>
  );
}