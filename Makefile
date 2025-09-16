CXX := g++
CXXFLAGS := -Wall -Wextra -std=c++17
INCLUDES = -Iclasses
BUILD_DIR := build

TASK_DIRS := $(wildcard 1/[0-9]*.*/)
TASKS := $(patsubst 1/%/,%,$(TASK_DIRS))

GENERAL_TASKS := $(filter-out 1.3,$(TASKS))

ALL_TARGETS := $(TASKS) 1.3.2

all: $(ALL_TARGETS)
	@echo "Все задания собраны!"

$(GENERAL_TASKS): %:
	@mkdir -p $(BUILD_DIR)
	@echo "Сборка задания $@..."
	@if [ -f "1/$@/$@.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/$@ 1/$@/$@.cpp; \
		echo "Собрано: $(BUILD_DIR)/$@ из 1/$@/$@.cpp"; \
	elif [ -f "1/$@/$(subst .,_,$@).cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/$@ 1/$@/$(subst .,_,$@).cpp; \
		echo "Собрано: $(BUILD_DIR)/$@ из 1/$@/$(subst .,_,$@).cpp"; \
	else \
		echo "Поиск любого .cpp файла в 1/$@..."; \
		CPP_FILE=$$(find "1/$@" -maxdepth 1 -name "*.cpp" | head -1); \
		if [ -n "$$CPP_FILE" ]; then \
			$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/$@ "$$CPP_FILE"; \
			echo "Собрано: $(BUILD_DIR)/$@ из $$CPP_FILE"; \
		else \
			echo "Ошибка: Не найден .cpp файл для задания $@"; \
		fi; \
	fi

1.3:
	@mkdir -p $(BUILD_DIR)
	@echo "Сборка задания 1.3 (основной)..."
	@if [ -f "1/1.3/1_3.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/1.3 1/1.3/1_3.cpp; \
		echo "Собрано: $(BUILD_DIR)/1.3 из 1/1.3/1_3.cpp"; \
	else \
		echo "Ошибка: Не найден файл 1/1.3/1_3.cpp"; \
	fi

1.3.2:
	@mkdir -p $(BUILD_DIR)
	@echo "Сборка подзадания 1.3.2..."
	@if [ -f "1/1.3/1_3_2.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/1.3.2 1/1.3/1_3_2.cpp; \
		echo "Собрано: $(BUILD_DIR)/1.3.2 из 1/1.3/1_3_2.cpp"; \
	else \
		echo "Ошибка: Не найден файл 1/1.3/1_3_2.cpp"; \
	fi

list:
	@echo "Доступные задания: $(TASKS)"
	@echo "Доступные подзадания: 1.3.2"
	@echo "Используйте: make <номер_задания>"
	@echo "Например: make 1.5 или make 1.3.2"
	@echo "Или: make all для сборки всего"

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all $(GENERAL_TASKS) 1.3 1.3.2 list clean