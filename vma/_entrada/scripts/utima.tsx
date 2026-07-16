import React, { useState } from 'react';
import { Download, ShieldCheck, TreeDeciduous, Lock, RefreshCw, Briefcase, CheckCircle, CreditCard, ScrollText, AlertTriangle } from 'lucide-react';

const DOWNLOAD_URL = "https://www.albmanvic.wordpress.com";

export default function K3Hub() {
  const [activeTab, setActiveTab] = useState('home');
  const [isVerifying, setIsVerifying] = useState(false);
  const [hashValue, setHashValue] = useState('');
  const [showVerification, setShowVerification] = useState(false);
  
  const [formData, setFormData] = useState({
    companyName: '',
    subject: '',
    taxId: ''
  });

  const contractText = `CONTRATO DE LICENCIA DE USO Y EXPLOTACIÓN VITALICIA
  
  LICENCIANTE: Víctor (Amor)
  LICENCIATARIO: ${formData.companyName || '[EMPRESA]'} (${formData.taxId || '[NIF/CIF]'})

  1. OBJETO: El Licenciante otorga una licencia no exclusiva, intransferible y de carácter vitalicio para el uso y explotación comercial de la tecnología "HashCode".
  2. PRECIO: El Licenciatario abona la cantidad única de 133,00 € mediante transferencia Bizum, previo proceso de verificación de identidad (0,33€ + 0,01€).
  3. AVAL LEGAL: La suma de las tres transacciones (0,33€, 0,01€ y 133,00€) capturadas por el Licenciatario constituye el sello de garantía y aval de la firma digital de Víctor (Amor).
  4. PROPIEDAD INTELECTUAL: La propiedad intelectual permanece bajo la titularidad del Licenciante.
  5. EXCLUSIÓN DE GARANTÍAS: El software se suministra "tal cual" (AS-IS).
  
  Firmado electrónicamente por: AMOR`;

  const handleVerify = () => {
    if (!formData.companyName || !formData.subject || !formData.taxId) return;

    setIsVerifying(true);
    setTimeout(() => {
      const inputSeed = formData.companyName + formData.taxId;
      const mockHash = "K3-SEC-" + btoa(inputSeed).substring(0, 8).toUpperCase() + "-" + 
                       Math.random().toString(16).substring(2, 6).toUpperCase();
                       
      setHashValue(mockHash);
      setIsVerifying(false);
      setShowVerification(true);
    }, 1500);
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
            <p className="text-xl text-slate-600 max-w-lg mx-auto">Eficiencia computacional. Acceso libre para investigación personal.</p>
            <a href={DOWNLOAD_URL} className="inline-flex items-center gap-3 bg-stone-900 text-white px-8 py-4 rounded-xl font-bold hover:bg-emerald-800 transition transform hover:scale-105 shadow-lg">
              <Download className="w-5 h-5" /> Descargar Software (.zip)
            </a>
          </div>
        ) : (
          <div className="space-y-8 animate-in fade-in slide-in-from-bottom-4 duration-500">
            <div className="bg-amber-50 border border-amber-200 p-4 rounded-lg flex gap-3 text-amber-800 text-sm">
              <AlertTriangle className="w-5 h-5 shrink-0" />
              <p><strong>Uso Exclusivo Empresas:</strong> Esta vía de formalización es solo para entidades con actividad comercial.</p>
            </div>

            {!showVerification ? (
              <div className="bg-white p-8 rounded-2xl border border-stone-200 shadow-sm space-y-6">
                <h4 className="font-bold text-lg text-emerald-900 flex items-center gap-2">
                  <ShieldCheck className="w-5 h-5" /> Iniciar Formalización de Licencia
                </h4>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <input type="text" placeholder="Nombre de la Empresa" className="w-full p-3 rounded-lg border border-stone-200" value={formData.companyName} onChange={(e) => setFormData({...formData, companyName: e.target.value})} />
                  <input type="text" placeholder="NIF/CIF" className="w-full p-3 rounded-lg border border-stone-200" value={formData.taxId} onChange={(e) => setFormData({...formData, taxId: e.target.value})} />
                </div>
                <input type="text" placeholder="Motivo de uso comercial" className="w-full p-3 rounded-lg border border-stone-200" value={formData.subject} onChange={(e) => setFormData({...formData, subject: e.target.value})} />
                <button onClick={handleVerify} disabled={isVerifying} className="w-full flex items-center justify-center gap-2 bg-emerald-600 text-white px-6 py-3 rounded-lg font-bold hover:bg-emerald-700 transition">
                  {isVerifying ? <RefreshCw className="animate-spin" /> : <Lock />} Generar Handshake
                </button>
              </div>
            ) : (
              <div className="space-y-6 animate-in fade-in">
                <div className="bg-white p-8 rounded-2xl border-2 border-emerald-500 shadow-xl space-y-4">
                    <h3 className="font-bold text-xl text-emerald-900 flex items-center gap-2"><ScrollText/> Contrato (Vigente tras pago)</h3>
                    <div className="font-mono text-xs bg-stone-50 p-4 rounded border border-stone-200 h-48 overflow-y-auto whitespace-pre-line text-slate-700">
                        {contractText}
                    </div>
                </div>

                <div className="bg-emerald-900 p-6 rounded-2xl text-emerald-100 space-y-4 shadow-xl">
                  <h3 className="font-bold text-xl text-white">Handshake Activo: {hashValue}</h3>
                  <div className="space-y-3 text-sm opacity-90 border-t border-emerald-700 pt-4">
                    <p><strong>Paso 1:</strong> Bizum <strong>0,33 €</strong> al +674 941 533. (Concepto: <strong>{hashValue}</strong>).</p>
                    <p><strong>Paso 2:</strong> Bizum <strong>0,01 €</strong> al +674 941 533. (Concepto: <strong>{hashValue}-VAL</strong>).</p>
                    <p><strong>Paso 3:</strong> Bizum <strong>133,00 €</strong> al +674 941 533. (Concepto: <strong>{hashValue}-LIC</strong>).</p>
                    <div className="bg-emerald-800 p-3 rounded border border-emerald-600 text-emerald-50 mt-4">
                       <CheckCircle className="inline w-4 h-4 mr-2"/> <strong>IMPORTANTE:</strong> La suma de estas 3 transacciones capturadas es su <strong>Aval Legal y Firma Digital de AMOR</strong>. Guárdelas como su licencia de por vida.
                    </div>
                  </div>
                </div>
              </div>
            )}
          </div>
        )}
      </main>
    </div>
  );
}