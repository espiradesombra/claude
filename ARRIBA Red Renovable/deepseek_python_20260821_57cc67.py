class MegabateriaConPerdidasNegativas:
    def __init__(self):
        # 400 Kilómetros con diferentes fases de patada
        self.kilometros = []
        for i in range(400):
            if i < 120:
                fases = 3  # 30%
            elif i < 280:
                fases = 5  # 40%
            else:
                fases = 7  # 30%
            self.kilometros.append(KilometroConPerdidasNegativas(i, 10, 15, 0.85, 0.90, 1.5, fases))
        
        # 150 módulos de molinos (100 con Quijote + 50 sin Quijote)
        self.modulos_molinos = []
        for i in range(150):
            if i < 100:
                self.modulos_molinos.append(ModuloMolinosConQuijote(i))
            else:
                self.modulos_molinos.append(ModuloMolinosSinQuijote(i))
        
        # Estado global
        self.energia_total_generada = 0.0
        self.energia_total_consumida = 0.0
        self.potencia_instantanea = 0.0
        
    def ejecutar_ciclo_completo(self, duracion=3600):
        """Ejecuta un ciclo completo de 1 hora."""
        # 1. Molinos generan energía (según viento variable)
        energia_molinos = 0
        for modulo in self.modulos_molinos:
            # Viento variable (simulación realista)
            viento = 12 + 4 * np.sin(np.random.random() * 2 * np.pi)
            energia_molinos += modulo.generar_energia(viento, duracion)
        
        # 2. Kilómetros generan energía (con pérdidas negativas)
        energia_KM = 0
        for km in self.kilometros:
            # Cada Kilómetro hace ciclos durante la hora
            ciclos = int(duracion / km.tiempo_ciclo)
            for _ in range(ciclos):
                E_balance, E_gen, E_cons = km.ciclo_con_patada(km.fases)
                energia_KM += E_balance
                self.energia_total_generada += E_gen
                self.energia_total_consumida += E_cons
        
        # 3. Balance total
        energia_total = energia_molinos + energia_KM
        self.potencia_instantanea = energia_total / duracion
        
        return {
            'energia_molinos': energia_molinos,
            'energia_KM': energia_KM,
            'energia_total': energia_total,
            'potencia_media': self.potencia_instantanea,
            'ciclos_KM': sum(km.ciclos for km in self.kilometros),
            'eficiencia_aparente': (energia_KM + energia_molinos) / energia_KM
        }