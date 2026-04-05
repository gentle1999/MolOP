# =============================================================================
# 🛠️ MolOP Makefile
# =============================================================================

PACKAGE_NAME := molop
REPOSITORY_URL := https://github.com/gentle1999/MolOP

# --- 自动检测版本号 ---
# 尝试通过 importlib 读取已安装包的版本，如果失败则显示 "dynamic"
VERSION := $(shell uv run python -c "from importlib.metadata import version; print(version('$(PACKAGE_NAME)'))" 2>/dev/null || echo "dynamic")

# 检测操作系统，用于打开浏览器命令
DETECTED_OS := $(shell uname)
ifeq ($(DETECTED_OS), Darwin)
	OPEN_CMD := open
else
	OPEN_CMD := xdg-open
endif

.PHONY: help init install-uv install sync update tree format format-check lint lint-check type-check pyright check-types check test test-cov clean distclean build release docker-build docker-up docker-down gen-typing-stubs check-typing-stubs gen-format-transform-stubs check-format-transform-stubs gen-cli-transform-stubs check-cli-transform-stubs docs-serve docs-build docs-build-strict

# =============================================================================
# 📝 帮助文档
# =============================================================================
help:
	@echo "📚 \033[1;34mMolOP Makefile Helper\033[0m"
	@echo ""
	@echo "📦 \033[1;33mDependency Management:\033[0m"
	@echo "  make install-uv  ⬇️ Install uv (The package manager)"
	@echo "  make init        🚀 Initialize environment (Interactive: Pin Python + Update Metadata)"
	@echo "  make update      🔄 Update dependencies (uv lock --upgrade)"
	@echo "  make tree        🌳 Show dependency tree"
	@echo ""
	@echo "🎨 \033[1;33mCode Quality:\033[0m"
	@echo "  make format      ✨ Format code (ruff format)"
	@echo "  make format-check ✅ Check formatting without changing files"
	@echo "  make lint        🔧 Lint code and apply safe fixes (ruff check --fix)"
	@echo "  make lint-check  🔍 Lint code without changing files"
	@echo "  make type-check  🦆 Static type check (mypy)"
	@echo "  make gen-typing-stubs     🧩 Generate typing stubs (.pyi)"
	@echo "  make check-typing-stubs   ✅ Verify typing stubs are up to date"
	@echo "  make gen-format-transform-stubs     🧩 Generate format transform stubs"
	@echo "  make check-format-transform-stubs   ✅ Verify format transform stubs are up to date"
	@echo "  make gen-cli-transform-stubs        🧩 Generate CLI transform stubs"
	@echo "  make check-cli-transform-stubs      ✅ Verify CLI transform stubs are up to date"
	@echo "  make check       🛡️ Run non-mutating quality checks"
	@echo ""
	@echo "🧪 \033[1;33mTesting:\033[0m"
	@echo "  make test        🌡️ Run unit tests"
	@echo "  make test-cov    📊 Run tests with HTML coverage report & open it"
	@echo ""
	@echo "🏗️ \033[1;33mBuild & Release:\033[0m"
	@echo "  make build       📦 Build package (sdist + wheel)"
	@echo "  make release     🚀 Print release instructions (Git Tag / GitHub Release)"
	@echo ""
	@echo "🧹 \033[1;33mUtilities:\033[0m"
	@echo "  make clean       🧹 Clean build artifacts & cache"
	@echo "  make distclean   🗑️ Clean EVERYTHING (including .venv)"
	@echo "  make rename      🏷️ [Template Tool] Global rename of '$(PACKAGE_NAME)'"
	@echo ""
	@echo "🐳 \033[1;33mDocker Integration:\033[0m"
	@echo "  make docker-build  🏗️ Build Docker image"
	@echo "  make docker-up     🚀 Run with Docker Compose"
	@echo "  make docker-down   🛑 Stop Docker services"
	@echo ""
	@echo "📖 \033[1;33mDocumentation:\033[0m"
	@echo "  make docs-serve        🌐 Serve documentation locally"
	@echo "  make docs-build        🏗️ Build documentation"
	@echo "  make docs-build-strict ✅ Build documentation with strict mode"
	@echo ""
	@echo "📌 Current Package: $(PACKAGE_NAME)"
	@echo "📌 Detected Version: $(VERSION)"

# =============================================================================
# 📦 依赖管理
# =============================================================================

# --- 🛠️ 安装 uv 工具 ---
install-uv:
	@echo "⬇️ Installing uv via official script..."
	@curl -LsSf https://astral.sh/uv/install.sh | sh
	@echo "✅ uv installed! You might need to restart your shell or run 'source $$HOME/.cargo/env'"

# --- 🚀 初始化项目 ---
init:
	@# 1. 检查 uv 是否存在
	@if ! command -v uv >/dev/null 2>&1; then \
		echo "⚠️ uv not found. Installing now..."; \
		$(MAKE) install-uv; \
	fi
	
	@# 2. 交互式询问 Python 版本
	@echo "🐍 \033[1;33mLet's configure your Python environment.\033[0m"
	@read -p "👉 Enter Python version to use (default: 3.10): " py_ver; \
	if [ -z "$$py_ver" ]; then \
		py_ver="3.10"; \
	fi; \
	\
	echo "📌 Pinning Python version to $$py_ver (.python-version)..."; \
	uv python pin $$py_ver; \
	\
	echo "📝 Updating pyproject.toml (requires-python >= $$py_ver)..."; \
	sed -i.bak "s/^requires-python = .*/requires-python = \">=$$py_ver\"/" pyproject.toml && rm pyproject.toml.bak

	@# 3. 同步依赖
	@echo "🚀 Installing dependencies..."
	uv sync --all-groups
	
	@echo "✅ Environment ready! Activate with: source .venv/bin/activate"

# 别名
install: init

# 仅仅同步
sync:
	uv sync

# 升级依赖
update:
	@echo "🔄 Updating dependencies..."
	uv lock --upgrade
	uv sync

# 显示依赖树
tree:
	uv tree

# =============================================================================
# 🎨 代码质量
# =============================================================================
format:
	@echo "🎨 Running Ruff Formatter..."
	uv run ruff format .

format-check:
	@echo "🎨 Checking Ruff formatting..."
	uv run ruff format . --check

lint:
	@echo "🔧 Running Ruff Linter with fixes..."
	uv run ruff check . --fix

lint-check:
	@echo "🔍 Running Ruff Linter..."
	uv run ruff check .

type-check:
	@echo "🦆 Running Mypy Type Checker..."
	uv run mypy src

pyright:
	@echo "🦆 Running Pyright Type Checker..."
	uv run pyright

check-types: type-check pyright

check: format-check lint-check check-types check-typing-stubs check-format-transform-stubs check-cli-transform-stubs

# =============================================================================
# 🧩 Typing stub generation
# =============================================================================

gen-typing-stubs:
	uv run python scripts/generate_io_typing_catalog.py

check-typing-stubs:
	uv run python scripts/generate_io_typing_catalog.py --check

gen-format-transform-stubs:
	uv run python scripts/generate_chemfile_format_transform_stubs.py

check-format-transform-stubs:
	uv run python scripts/generate_chemfile_format_transform_stubs.py --check

gen-cli-transform-stubs:
	uv run python scripts/generate_cli_transform_stubs.py

check-cli-transform-stubs:
	uv run python scripts/generate_cli_transform_stubs.py --check

# =============================================================================
# 🧪 测试与覆盖率
# =============================================================================
test:
	@echo "🧪 Running Pytest..."
	uv run pytest

test-cov:
	@echo "📊 Running Test Coverage..."
	uv run pytest --cov=src --cov-report=html
	@echo "🌍 Opening coverage report..."
	@$(OPEN_CMD) htmlcov/index.html

# =============================================================================
# 🏗️ 构建与打包
# =============================================================================
build: clean check test
	@echo "🏗️ Building package (Hatchling + UV)..."
	uv build

# =============================================================================
# 🚀 发布流程
# =============================================================================
release:
	@echo ""
	@echo "🚀 \033[1;32mReady to release version?\033[0m (Current detection: $(VERSION))"
	@echo "---------------------------------------------------"
	@echo "Since you are using CI/CD driven releases:"
	@echo "1. Commit all your changes."
	@echo "2. Go to your repository release page:"
	@echo "   🔗 GitHub: $(REPOSITORY_URL)/releases/new"
	@echo "3. Draft a new release and create a tag (e.g., v0.1.0)."
	@echo "---------------------------------------------------"

# =============================================================================
# 🧹 清理
# =============================================================================
clean:
	@echo "🧹 Cleaning artifacts..."
	rm -rf dist build htmlcov coverage.xml .coverage
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name "*.egg-info" -exec rm -rf {} +
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
	find . -type d -name ".ruff_cache" -exec rm -rf {} +
	find . -type d -name ".mypy_cache" -exec rm -rf {} +

distclean: clean
	@echo "🗑️ Removing virtual environment (.venv)..."
	rm -rf .venv .python-version
	@echo "✨ Project is clean. Run 'make init' to restart."

# =============================================================================
# 🐳 Docker 常用命令
# =============================================================================
docker-build:
	@echo "🏗️ Building Docker image for $(PACKAGE_NAME)..."
	docker build -t $(PACKAGE_NAME):latest .

docker-up:
	@echo "🚀 Starting services..."
	docker compose up -d --build
	@echo "📜 Use 'docker compose logs -f' to follow logs."

docker-down:
	@echo "🛑 Stopping services..."
	docker compose down

# =============================================================================
# 📖 文档相关
# =============================================================================
docs-serve:
	@echo "🌐 Serving documentation..."
	uv run mkdocs serve -f mkdocs.dev.yml -a 127.0.0.1:8000

docs-build:
	@echo "🏗️ Building documentation..."
	uv run mkdocs build

docs-build-strict:
	@echo "✅ Building documentation (strict mode)..."
	uv run mkdocs build --strict
