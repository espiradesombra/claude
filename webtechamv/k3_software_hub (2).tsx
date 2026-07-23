import React, { useState } from 'react';
import { Download, ShieldCheck, TreeDeciduous, Lock, RefreshCw, Briefcase, FileText, AlertTriangle, CheckCircle } from 'lucide-react';

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

  const handleVerify = () => {
    if (!formData.companyName || !formData.subject || !formData.taxId) {
      alert("Por favor, rellene todos los campos: Nombre de Empresa, NIF/CIF y Motivo de uso.");
      return;
    }

    setIsVerifying(true);
    // Simulación de proceso de Handshake
    setTimeout(() => {
      const inputSeed = formData.companyName + formData.taxId;
      const mockHash = "K3-SEC-" + btoa(inputSeed).substring(0, 16).toUpperCase() + "-" + 
                       Math.random().toString(16).substring(2, 6).toUpperCase();
                       
      setHashValue(mockHash);
      setIsVerifying(false);
      setShowVerification(true);
    }, 2000);
  };

  return (
    <div className="min-h-screen bg-stone-50 text-slate-900 font-sans">
      {/* Navigation */}
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
              Eficiencia computacional. Acceso libre para investigación personal.
            </p>
            <a 
              href={DOWNLOAD_URL}
              className="inline-flex items-center gap-3 bg-stone-900 text-white px-8 py-4 rounded-xl font-bold hover:bg-emerald-800 transition transform hover:scale-105 shadow-lg"
            >
              <Download className="w-5 h-5" /> Descargar Software (.zip)
            </a>
          </div>
        ) : (
          <div className="space-y-8 animate-in fade-in slide-in-from-bottom-4 duration-500">
            <div className="bg-amber-50 border border-amber-200 p-4 rounded-lg flex gap-3 text-amber-800 text-sm">
              <AlertTriangle className="w-5 h-5 shrink-0" />
              <p><strong>Aviso Importante:</strong> Esta sección es exclusivamente para uso comercial. No proceda si usted es un usuario particular, utilice la descarga gratuita en la sección "Acceso".</p>
            </div>

            <header>
              <h2 className="text-3xl font-bold text-emerald-900 mb-2">Protocolo de Licencia Vitalicia</h2>
              <p className="text-slate-600">Registro y validación de contrato para entidades comerciales.</p>
            </header>

            {!showVerification ? (
              <div className="bg-white p-8 rounded-2xl border border-stone-200 shadow-sm space-y-6">
                <h4 className="font-bold text-lg text-emerald-900 flex items-center gap-2">
                  <ShieldCheck className="w-5 h-5" /> Iniciar Verificación Comercial
                </h4>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <input type="text" placeholder="Nombre de la Empresa" className="w-full p-3 rounded-lg border border-stone-200" value={formData.companyName} onChange={(e) => setFormData({...formData, companyName: e.target.value})} />
                  <input type="text" placeholder="NIF/CIF" className="w-full p-3 rounded-lg border border-stone-200" value={formData.taxId} onChange={(e) => setFormData({...formData, taxId: e.target.value})} />
                </div>
                <input type="text" placeholder="Motivo de uso comercial (Ej: Integración en Antivirus)" className="w-full p-3 rounded-lg border border-stone-200" value={formData.subject} onChange={(e) => setFormData({...formData, subject: e.target.value})} />
                <button onClick={handleVerify} disabled={isVerifying} className="w-full flex items-center justify-center gap-2 bg-emerald-600 text-white px-6 py-3 rounded-lg font-bold hover:bg-emerald-700 transition">
                  {isVerifying ? <RefreshCw className="animate-spin" /> : <Lock />} Generar Handshake de Validación
                </button>
              </div>
            ) : (
              <div className="space-y-6 animate-in fade-in">
                <div className="bg-emerald-900 p-6 rounded-2xl text-emerald-100 space-y-4 shadow-xl">
                  <h3 className="font-bold text-xl text-white">Estado: Validado técnicamente</h3>
                  <div className="font-mono text-sm bg-black/30 p-4 rounded-lg break-all border border-emerald-700">
                    {hashValue}
                  </div>
                  <p className="text-sm opacity-90">Entidad: <strong>{formData.companyName}</strong> ({formData.taxId})</p>
                </div>

                <div className="bg-white p-6 rounded-2xl border border-emerald-200 shadow-sm space-y-4">
                    <h4 className="font-bold text-emerald-900 flex items-center gap-2"><CheckCircle className="text-emerald-500"/> Instrucciones de Formalización</h4>
                    <p className="text-sm text-slate-700">Para finalizar la emisión del contrato firmado por <strong>Amor</strong>:</p>
                    <ol className="list-decimal pl-5 space-y-2 text-sm text-slate-700">
                        <li>Realice transferencia de <strong>0,33 €</strong> vía Bizum al <strong>+674 941 533</strong>. Concepto: <strong>{hashValue.substring(0, 8)}</strong>.</li>
                        <li>Realice una devolución de <strong>0,01 €</strong> al mismo número. Esto confirma su titularidad bancaria.</li>
                        <li>Tras la recepción, le enviaremos por email el Contrato de Licencia Vitalicia completo para su firma digital.</li>
                    </ol>
                </div>

                <div className="bg-stone-100 p-6 rounded-xl border border-stone-200 text-xs text-slate-500 font-mono">
                    <h5 className="font-bold mb-2">CONTRATO DE LICENCIA (Extracto Legal)</h5>
                    <p>LICENCIANTE: Víctor (Amor) | LICENCIATARIO: {formData.companyName} ({formData.taxId})</p>
                    <p className="mt-2">Objeto: Licencia vitalicia de uso y explotación sobre algoritmos de identificación HashCode.</p>
                    <p className="mt-2">Condiciones: El software se entrega "AS-IS". La validez de este contrato está sujeta a la recepción del pago de validación bancaria y la aceptación de los términos de uso comercial.</p>
                </div>
              </div>
            )}
          </div>
        )}
      </main>
    </div>
  );
}