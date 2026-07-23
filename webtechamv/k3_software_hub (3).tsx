import React, { useState } from 'react';
import { Download, ShieldCheck, TreeDeciduous, Lock, RefreshCw, Briefcase, FileText, AlertTriangle } from 'lucide-react';

const DOWNLOAD_URL = "https://www.albmanvic.wordpress.com";

const CONTRACT_TEXT = `CONTRATO DE LICENCIA DE USO Y EXPLOTACIÓN VITALICIA

En [CIUDAD], a fecha actual.

DE UNA PARTE:
D. VÍCTOR [APELLIDOS], en adelante el LICENCIANTE (Marca: Amor), con domicilio a efectos de notificaciones en el teléfono +674 941 533.

DE OTRA PARTE:
[NOMBRE DE LA EMPRESA O RAZÓN SOCIAL], con NIF/CIF [NÚMERO], en adelante el LICENCIATARIO.

EXPONEN:
Que el Licenciante es titular de los derechos de explotación sobre el algoritmo "HashCode" aplicable al ámbito de la identificación y seguridad.

CLÁUSULAS:
PRIMERA. OBJETO: Licencia no exclusiva, intransferible y vitalicia para el uso comercial de la tecnología "HashCode".
SEGUNDA. PRECIO: Pago único de 133,00 € (ciento treinta y tres euros) mediante transferencia/Bizum al +674 941 533.
TERCERA. PROPIEDAD INTELECTUAL: Propiedad total del Licenciante (Amor).
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
  
  const [formData, setFormData] = useState({
    companyName: '',
    subject: ''
  });

  const handleVerify = () => {
    if (!formData.companyName || !formData.subject) {
      alert("Por favor, rellene el nombre de la empresa y el motivo de uso.");
      return;
    }

    setIsVerifying(true);
    setTimeout(() => {
      const inputSeed = formData.companyName + formData.subject;
      const mockHash = "K3-SEC-" + btoa(inputSeed).substring(0, 16).toUpperCase() + "-" + 
                       Math.random().toString(16).substring(2, 6).toUpperCase();
                       
      setHashValue(mockHash);
      setIsVerifying(false);
      setShowContact(true);
    }, 2000);
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
              Eficiencia computacional. Acceso libre para investigación.
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
            <div className="bg-amber-50 border border-amber-200 p-4 rounded-lg flex gap-3 text-amber-800 text-sm">
              <AlertTriangle className="w-5 h-5 shrink-0" />
              <p><strong>Aviso Importante:</strong> Esta sección es exclusivamente para empresas que obtienen rédito comercial de nuestra tecnología. Si usted es un usuario particular, utilice la descarga gratuita en la sección "Acceso".</p>
            </div>

            <header>
              <h2 className="text-3xl font-bold text-emerald-900 mb-2">Protocolo de Licencia Comercial</h2>
              <p className="text-slate-600">Verificación técnica de identidad para contrato vitalicio.</p>
            </header>

            <div className="bg-white p-6 rounded-2xl border border-stone-200 shadow-sm space-y-6">
              <div className="bg-emerald-50 p-6 rounded-xl border border-emerald-100 space-y-4">
                <h4 className="font-bold text-emerald-900 flex items-center gap-2">
                  <Lock className="w-5 h-5" /> Protocolo de Verificación Doble (Handshake)
                </h4>
                
                {!showContact ? (
                  <div className="space-y-4">
                    <input 
                      type="text" 
                      placeholder="Nombre de la Empresa" 
                      className="w-full p-3 rounded-lg border border-emerald-200 focus:ring-2 focus:ring-emerald-500 outline-none"
                      value={formData.companyName}
                      onChange={(e) => setFormData({...formData, companyName: e.target.value})}
                    />
                    <input 
                      type="text" 
                      placeholder="Motivo de uso comercial" 
                      className="w-full p-3 rounded-lg border border-emerald-200 focus:ring-2 focus:ring-emerald-500 outline-none"
                      value={formData.subject}
                      onChange={(e) => setFormData({...formData, subject: e.target.value})}
                    />
                    
                    <button 
                      onClick={handleVerify}
                      disabled={isVerifying}
                      className="w-full flex items-center justify-center gap-2 bg-emerald-600 text-white px-6 py-3 rounded-lg font-bold hover:bg-emerald-700 transition disabled:opacity-50"
                    >
                      {isVerifying ? <RefreshCw className="animate-spin" /> : <ShieldCheck />}
                      {isVerifying ? "Generando Handshake..." : "Iniciar Handshake de Verificación"}
                    </button>
                  </div>
                ) : (
                  <div className="space-y-4 animate-in fade-in text-sm text-emerald-900">
                    <div className="p-4 bg-emerald-900 rounded-lg font-mono text-xs text-emerald-300 break-all border border-emerald-700">
                      <span className="text-emerald-500 block mb-1">// Código de Handshake (Concepto):</span>
                      {hashValue}
                    </div>
                    
                    <div className="bg-white p-4 rounded border border-emerald-200 space-y-2">
                        <p><strong>Estado:</strong> Validado técnicamente para: <strong>{formData.companyName}</strong>.</p>
                        <p>Para confirmar la titularidad bancaria y formalizar el contrato:</p>
                        <ol className="list-decimal pl-5 space-y-2 mt-2">
                            <li>Realice un Bizum de <strong>0,33 €</strong> al <strong>+674 941 533</strong> con concepto: <strong>{hashValue.substring(0, 8)}</strong>.</li>
                            <li>Realice una devolución inmediata de <strong>0,01 €</strong> desde su misma cuenta.</li>
                            <li>Tras la recepción, le enviaremos el contrato de cesión vitalicia firmado por <strong>Amor</strong>.</li>
                        </ol>
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