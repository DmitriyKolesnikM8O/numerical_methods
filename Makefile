# Компилятор и флаги
CXX := g++
CXXFLAGS := -Wall -Wextra -std=c++17
INCLUDES = -Itasks/classes
BUILD_DIR := build
TASKS_ROOT := tasks

# Находим все подпапки в директории Tasks/1
TASK_DIRS := $(wildcard $(TASKS_ROOT)/1/[0-9]*.*/)
TASKS := $(patsubst $(TASKS_ROOT)/1/%/,%,$(TASK_DIRS))

# Исключаем 1.3 из общих задач, так как для него есть отдельное правило
GENERAL_TASKS := $(filter-out 1.3,$(TASKS))

# Все цели для сборки (основные задания + подзадание 1.3.2)
ALL_TARGETS := $(TASKS) 1.3.2

# Правило по умолчанию - собираем все
all: $(ALL_TARGETS)
	@echo "Все задания собраны!"

# Общее правило для большинства заданий
$(GENERAL_TASKS): %:
	@mkdir -p $(BUILD_DIR)
	@echo "Сборка задания $@..."
	@if [ -f "$(TASKS_ROOT)/1/$@/$@.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/$@ $(TASKS_ROOT)/1/$@/$@.cpp; \
		echo "Собрано: $(BUILD_DIR)/$@ из $(TASKS_ROOT)/1/$@/$@.cpp"; \
	elif [ -f "$(TASKS_ROOT)/1/$@/$(subst .,_,$@).cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/$@ $(TASKS_ROOT)/1/$@/$(subst .,_,$@).cpp; \
		echo "Собрано: $(BUILD_DIR)/$@ из $(TASKS_ROOT)/1/$@/$(subst .,_,$@).cpp"; \
	else \
		echo "Поиск любого .cpp файла в $(TASKS_ROOT)/1/$@..."; \
		CPP_FILE=$$(find "$(TASKS_ROOT)/1/$@" -maxdepth 1 -name "*.cpp" | head -1); \
		if [ -n "$$CPP_FILE" ]; then \
			$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/$@ "$$CPP_FILE"; \
			echo "Собрано: $(BUILD_DIR)/$@ из $$CPP_FILE"; \
		else \
			echo "Ошибка: Не найден .cpp файл для задания $@"; \
		fi; \
	fi

# Конкретное правило для 1.3 (основной файл)
1.3:
	@mkdir -p $(BUILD_DIR)
	@echo "Сборка задания 1.3 (основной)..."
	@if [ -f "$(TASKS_ROOT)/1/1.3/1_3.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/1.3 $(TASKS_ROOT)/1/1.3/1_3.cpp; \
		echo "Собрано: $(BUILD_DIR)/1.3 из $(TASKS_ROOT)/1/1.3/1_3.cpp"; \
	else \
		echo "Ошибка: Не найден файл $(TASKS_ROOT)/1/1.3/1_3.cpp"; \
	fi

# Правило для подзадания 1.3.2
1.3.2:
	@mkdir -p $(BUILD_DIR)
	@echo "Сборка подзадания 1.3.2..."
	@if [ -f "$(TASKS_ROOT)/1/1.3/1_3_2.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/1.3.2 $(TASKS_ROOT)/1/1.3/1_3_2.cpp; \
		echo "Собрано: $(BUILD_DIR)/1.3.2 из $(TASKS_ROOT)/1/1.3/1_3_2.cpp"; \
	else \
		echo "Ошибка: Не найден файл $(TASKS_ROOT)/1/1.3/1_3_2.cpp"; \
	fi

# Показать доступные цели
list:
	@echo "Доступные задания: $(TASKS)"
	@echo "Доступные подзадания: 1.3.2"
	@echo "Используйте: make <номер_задания>"
	@echo "Например: make 1.5 или make 1.3.2"
	@echo "Или: make all для сборки всего"

# Очистка
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all $(GENERAL_TASKS) 1.3 1.3.2 list clean