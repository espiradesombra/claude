"""Plugin Manager — registro, métricas, discovery (UPS v0.0.065)."""

from __future__ import annotations

import importlib
import pkgutil
from dataclasses import dataclass, field
from time import perf_counter

from .plugin import IPlugin, PluginInfo


@dataclass
class PluginStats:
    executions: int = 0
    total_time_s: float = 0.0
    reuse_hits: int = 0
    failures: int = 0

    @property
    def average_time_ms(self) -> float:
        return (self.total_time_s / self.executions * 1000) if self.executions else 0.0


@dataclass
class PluginManager:
    plugins: dict[str, IPlugin] = field(default_factory=dict)
    stats: dict[str, PluginStats] = field(default_factory=dict)
    api_version: int = 1

    def register(self, plugin: IPlugin) -> None:
        info = plugin.info()
        if info.api_version > self.api_version:
            raise ValueError(f"plugin {info.name} API {info.api_version} > runtime {self.api_version}")
        self.plugins[info.name] = plugin
        self.stats.setdefault(info.name, PluginStats())

    def get(self, name: str) -> IPlugin:
        return self.plugins[name]

    def execute_timed(self, name: str, ctx) -> tuple[object, float]:
        plugin = self.get(name)
        st = self.stats[name]
        t0 = perf_counter()
        try:
            result = plugin.execute(ctx)
            st.executions += 1
            st.total_time_s += perf_counter() - t0
            return result, perf_counter() - t0
        except Exception:
            st.failures += 1
            raise

    def record_reuse(self, name: str) -> None:
        self.stats[name].reuse_hits += 1

    def discover_builtin(self) -> int:
        """Carga plugins del paquete plugins/."""
        import plugins as plugins_pkg

        count = 0
        for mod in pkgutil.iter_modules(plugins_pkg.__path__, plugins_pkg.__name__ + "."):
            if mod.name.endswith("_plugin"):
                module = importlib.import_module(mod.name)
                for attr in dir(module):
                    obj = getattr(module, attr)
                    if isinstance(obj, type) and issubclass(obj, IPlugin) and obj is not IPlugin:
                        try:
                            self.register(obj())
                            count += 1
                            break
                        except TypeError:
                            pass
        return count