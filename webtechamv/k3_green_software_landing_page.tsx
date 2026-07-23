import React, { useState } from 'react';
import { Leaf, ShieldCheck, Download, CreditCard, Code, TreeDeciduous } from 'lucide-react';

export default function App() {
  const [termsAccepted, setTermsAccepted] = useState(false);

  return (
    <div className="min-h-screen bg-stone-50 text-slate-900 font-sans selection:bg-emerald-200">
      
      {/* Navbar */}
      <nav className="p-6 flex justify-between items-center max-w-6xl mx-auto">
        <div className="flex items-center gap-2">
          <TreeDeciduous className="text-emerald-700 h-8 w-8" />
          <span className="text-2xl font-bold tracking-tight text-emerald-900">K3 Software</span>
        </div>
        <a href="#licensing" className="text-sm font-semibold text-emerald-800 hover:text-emerald-600 transition">
          Download & License
        </a>
      </nav>

      {/* Hero Section */}
      <header className="py-20 px-6 text-center max-w-4xl mx-auto">
        <div className="inline-flex items-center gap-2 px-4 py-1.5 bg-emerald-100 text-emerald-800 rounded-full text-sm font-medium mb-6">
          <Leaf className="w-4 h-4" /> Green Computing / Efficient Algorithms
        </div>
        <h1 className="text-5xl md:text-6xl font-extrabold text-slate-900 mb-6 leading-tight">
          Algoritmos <span className="text-emerald-700">inteligentes</span> para un futuro sostenible.
        </h1>
        <p className="text-lg text-slate-600 mb-10 max-w-2xl mx-auto">
          K3 reduce la carga computacional, optimizando recursos y disminuyendo la huella de carbono de tus sistemas. Eficiencia real para un mundo real.
        </p>
      </header>

      {/* Features Grid */}
      <section className="py-16 px-6 max-w-6xl mx-auto grid md:grid-cols-3 gap-8">
        {[
          { icon: Code, title: "Optimizado", desc: "Reducción drástica de ciclos de CPU." },
          { icon: ShieldCheck, title: "Seguridad Robusta", desc: "Cifrado y hashing de grado empresarial." },
          { icon: TreeDeciduous, title: "Ecodiseño", desc: "Menor consumo energético en cada proceso." }
        ].map((feat, i) => (
          <div key={i} className="bg-white p-8 rounded-2xl border border-stone-200 shadow-sm hover:shadow-md transition">
            <feat.icon className="w-10 h-10 text-emerald-600 mb-4" />
            <h3 className="text-xl font-bold mb-2">{feat.title}</h3>
            <p className="text-slate-600">{feat.desc}</p>
          </div>
        ))}
      </section>

      {/* Licensing / Download Section */}
      <section id="licensing" className="py-20 px-6 bg-emerald-900 text-white rounded-t-[3rem]">
        <div className="max-w-4xl mx-auto">
          <h2 className="text-4xl font-bold mb-12 text-center">Licencias y Descarga</h2>
          
          <div className="grid md:grid-cols-2 gap-8">
            {/* Free Tier */}
            <div className="bg-emerald-800/50 p-8 rounded-2xl border border-emerald-700">
              <h3 className="text-2xl font-bold mb-4">Uso Personal</h3>
              <p className="text-emerald-100 mb-6">Ideal para estudiantes, investigadores y proyectos open-source. Hashea y prueba sin costo.</p>
              <button className="w-full flex items-center justify-center gap-2 bg-emerald-600 hover:bg-emerald-500 py-3 rounded-lg font-bold transition">
                <Download className="w-5 h-5" /> Descarga Gratis (.zip)
              </button>
            </div>

            {/* Commercial Tier */}
            <div className="bg-white p-8 rounded-2xl text-slate-900 relative shadow-2xl">
              <span className="absolute -top-3 -right-3 bg-amber-400 text-amber-900 px-4 py-1 rounded-full text-sm font-bold">RECOMENDADO</span>
              <h3 className="text-2xl font-bold mb-4">Licencia Comercial</h3>
              <div className="text-4xl font-extrabold mb-6">133€ <span className="text-base font-normal text-slate-500">pago único</span></div>
              <p className="text-slate-600 mb-6">Permiso para uso en entornos productivos, integración en antivirus y software de terceros.</p>
              
              {/* Terms Checkbox */}
              <label className="flex items-start gap-3 mb-6 cursor-pointer">
                <input 
                  type="checkbox" 
                  className="mt-1 accent-emerald-600"
                  checked={termsAccepted}
                  onChange={(e) => setTermsAccepted(e.target.checked)}
                />
                <span className="text-sm text-slate-600">
                  He leído y acepto los <a href="#" className="underline font-bold">términos de la Licencia K3</a>. El pago constituye la firma digital del contrato.
                </span>
              </label>

              <button 
                disabled={!termsAccepted}
                className={`w-full flex items-center justify-center gap-2 py-3 rounded-lg font-bold transition ${
                  termsAccepted 
                    ? 'bg-emerald-700 text-white hover:bg-emerald-800' 
                    : 'bg-stone-200 text-stone-500 cursor-not-allowed'
                }`}
              >
                <CreditCard className="w-5 h-5" /> Adquirir Licencia (133€)
              </button>
            </div>
          </div>
        </div>
      </section>

      {/* Footer */}
      <footer className="py-8 text-center text-stone-500 text-sm">
        &copy; {new Date().getFullYear()} Víctor. Todos los derechos reservados.
      </footer>
    </div>
  );
}